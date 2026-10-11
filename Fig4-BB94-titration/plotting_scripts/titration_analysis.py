"""Shared gland-level statistics, reports and publication plots for Fig. 4."""
from __future__ import annotations

import argparse
import hashlib
import html
import importlib.metadata
import itertools
import json
from pathlib import Path
import platform
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def significance(p):
    if not np.isfinite(p) or not 0 <= p <= 1:
        raise ValueError("Significance requires a finite p-value between 0 and 1.")
    for threshold, label in [(1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (.05, "*")]:
        if p < threshold:
            return label
    return "ns"


def analyze(data, method="mann-whitney"):
    """Compare doses with DMSO by Mann-Whitney/Holm; legacy parametric modes use all pairs."""
    required = {"dose_um", "value", "cohort", "gland_id"}
    if not required <= set(data):
        raise ValueError(f"Missing columns: {sorted(required - set(data))}")
    if data[list(required)].isna().any().any():
        raise ValueError("Missing gland identity, dose or endpoint.")
    if data.duplicated(["cohort", "dose_um", "gland_id"]).any():
        raise ValueError("Each gland must contribute exactly one endpoint.")
    if not np.isfinite(data[["dose_um", "value"]].to_numpy(float)).all():
        raise ValueError("Doses and endpoints must be finite.")
    doses = sorted(data.dose_um.unique())
    if len(doses) < 2 or doses[0] != 0:
        raise ValueError("At least two doses, including vehicle (dose 0), are required.")
    groups = [data.loc[data.dose_um == dose, "value"].to_numpy(float) for dose in doses]
    if min(map(len, groups)) < 2:
        raise ValueError("Every dose needs at least two independent glands.")
    summary = data.groupby("dose_um").value.agg(n="count", mean="mean", sd="std").reset_index()
    summary["variance"] = summary.sd ** 2
    if method == "mann-whitney":
        rows = []
        for i, j in ((0, j) for j in range(1, len(groups))):
            a, b = groups[i], groups[j]
            pooled = np.concatenate([a, b])
            ties = len(np.unique(pooled)) != len(pooled)
            p_method = "exact" if min(len(a), len(b)) <= 8 and not ties else "asymptotic"
            test = stats.mannwhitneyu(b, a, alternative="two-sided", method=p_method, use_continuity=True)
            rows.append({"dose_a_um": doses[i], "dose_b_um": doses[j], "n_a": len(a), "n_b": len(b),
                         "mean_difference_b_minus_a": float(b.mean()-a.mean()),
                         "U_b": float(test.statistic),
                         "rank_biserial_b_minus_a": float(2*test.statistic/(len(a)*len(b))-1),
                         "p_raw": float(test.pvalue), "p_method": p_method, "ties_present": ties,
                         "method": "Two-sided Mann-Whitney U + Holm"})
        comparisons = pd.DataFrame(rows)
        comparisons["p_adjusted"] = multipletests(comparisons.p_raw, alpha=.05, method="holm")[1]
        comparisons["annotation"] = comparisons.p_adjusted.map(significance)
        diagnostics = {"method_selection": "Exact when either group has at most 8 glands and there are no pooled ties; otherwise asymptotic with tie and continuity correction.",
                       "comparison_family": "Each nonzero dose versus DMSO (0); Holm correction within this outcome.",
                       "p_method_counts": comparisons.p_method.value_counts().to_dict(),
                       "comparisons_with_ties": int(comparisons.ties_present.sum()),
                       "interpretation": "Rank-based independent-group comparison; not a test of mean difference or, without common distribution shape, a pure median test.",
                       "rank_effect": "2*U_b/(n_a*n_b)-1: positive means larger values tend to occur in dose B; ties contribute half to U."}
        return summary, comparisons, diagnostics
    if any(np.var(x, ddof=1) == 0 for x in groups):
        raise ValueError("A dose has zero variance; inspect inputs before parametric inference.")
    bf = stats.levene(*groups, center="median")
    diagnostics = {
        "brown_forsythe_statistic": float(bf.statistic),
        "brown_forsythe_p": float(bf.pvalue),
        "largest_to_smallest_variance": float(summary.variance.max() / summary.variance.min()),
        "small_sample_note": "Diagnostics do not establish normality or equal variances; small groups have limited power.",
        "shapiro_by_dose": [],
    }
    for dose, values in zip(doses, groups):
        if len(values) >= 3:
            sh = stats.shapiro(values)
            diagnostics["shapiro_by_dose"].append({"dose_um": float(dose), "n": len(values), "W": float(sh.statistic), "p": float(sh.pvalue)})
    omnibus = stats.f_oneway(*groups)
    diagnostics["ordinary_one_way_anova"] = {
        "F": float(omnibus.statistic), "p": float(omnibus.pvalue),
        "df_between": len(groups)-1, "df_within": len(data)-len(groups),
        "note": "Descriptive companion to Tukey; not a gate for the all-pairs procedure.",
    }
    rows = []
    if method == "tukey":
        test = stats.tukey_hsd(*groups)
        ci = test.confidence_interval(confidence_level=.95)
        for i, j in itertools.combinations(range(len(groups)), 2):
            rows.append({"dose_a_um": doses[i], "dose_b_um": doses[j],
                         "n_a": len(groups[i]), "n_b": len(groups[j]),
                         "difference_b_minus_a": float(test.statistic[j, i]),
                         "ci95_low": float(ci.low[j, i]), "ci95_high": float(ci.high[j, i]),
                         "p_adjusted": float(test.pvalue[j, i]), "method": "Tukey HSD / Tukey-Kramer"})
    elif method == "games-howell":
        # Games-Howell: studentized-range distribution with Welch-Satterthwaite df.
        k = len(groups)
        for i, j in itertools.combinations(range(k), 2):
            a, b = groups[i], groups[j]
            va, vb = np.var(a, ddof=1)/len(a), np.var(b, ddof=1)/len(b)
            se2 = va + vb
            df = se2**2 / (va**2/(len(a)-1) + vb**2/(len(b)-1))
            difference = float(b.mean() - a.mean())
            q = abs(difference) / np.sqrt(se2/2)
            p = float(stats.studentized_range.sf(q, k, df))
            half_ci = float(stats.studentized_range.ppf(.95, k, df) * np.sqrt(se2/2))
            rows.append({"dose_a_um": doses[i], "dose_b_um": doses[j], "n_a": len(a), "n_b": len(b),
                         "difference_b_minus_a": difference, "ci95_low": difference-half_ci,
                         "ci95_high": difference+half_ci, "p_adjusted": p,
                         "df": float(df), "method": "Games-Howell"})
    else:
        raise ValueError(f"Unknown method: {method}")
    comparisons = pd.DataFrame(rows)
    comparisons["annotation"] = comparisons.p_adjusted.map(significance)
    return summary, comparisons, diagnostics


def plot(data, summary, comparisons, kind, path, formats, font="Arial", width=1.9, height=1.75):
    import seaborn as sns
    font_manager.findfont(font_manager.FontProperties(family=font), fallback_to_default=False)
    doses = list(summary.dose_um)
    labels = [f"{d:g}" for d in doses]
    rc = {"font.family": font, "font.size": 7, "svg.fonttype": "none", "pdf.fonttype": 42,
          "axes.linewidth": .7, "xtick.major.width": .7, "ytick.major.width": .7}
    with matplotlib.rc_context(rc):
        fig, ax = plt.subplots(figsize=(width, height), dpi=300)
        sns.swarmplot(data=data, x="dose_um", y="value", order=doses, color="black",
                      size=1.9, linewidth=0, ax=ax)
        fig.canvas.draw()
        for dose, collection in zip(doses, ax.collections):
            collection.set_facecolor("black" if dose == 0 else "red")
        for i, row in enumerate(summary.itertuples()):
            ax.errorbar(i, row.mean, yerr=row.sd, fmt="_", color="gray", markersize=5,
                        elinewidth=.7, capsize=2, capthick=.7, zorder=5)
        highest = max(float(data.value.max()), float((summary["mean"] + summary.sd).max()))
        lowest = min(float(data.value.min()), float((summary["mean"] - summary.sd).min()))
        span = max(highest-lowest, 1)
        annotation_y = highest + .08*span
        ctrl = comparisons.loc[comparisons.dose_a_um == 0].set_index("dose_b_um")
        for i, dose in enumerate(doses[1:], start=1):
            annotation = ctrl.loc[dose, "annotation"]
            ax.text(i, annotation_y, annotation, ha="center", va="bottom",
                    fontsize=7 if annotation == "ns" else 8)
        if kind == "bud":
            baseline = int(data.baseline_day.iloc[0])
            final = int(data.final_day.iloc[0])
            ax.set_ylabel(r"$\mathregular{Log_2}$" + f"(bud ratio d{final}/d{baseline})")
            ax.set_ylim(min(-1, lowest-.1*span), max(4, annotation_y+.16*span))
            ax.set_yticks(np.arange(0, ax.get_ylim()[1]+.01, 2))
        else:
            ax.set_ylabel("Collagen IV intensity")
            ax.set_ylim(min(0, lowest-.1*span) if lowest < 0 else 0, max(4, annotation_y+.16*span))
            ax.set_yticks(np.arange(0, ax.get_ylim()[1]+.01, 1))
        ax.set_xticks(range(len(doses)), labels=labels, rotation=45, ha="right")
        ax.set_xlabel("BB-94 (µM)")
        ax.tick_params(axis="both", pad=2, length=3)
        ax.spines[["right", "top"]].set_visible(False)
        fig.tight_layout(pad=.25)
        for suffix in formats:
            fig.savefig(path.parent/(path.name+"."+suffix), dpi=450, bbox_inches="tight", pad_inches=.03)
        plt.close(fig)


def reports(folder, prefix, kind, data, summary, comparisons, diagnostics, audit, method, formats, alternative, cohort_sensitivity, cohort_metadata):
    title = "BB-94 titration: " + ("bud-count ratios" if kind == "bud" else "collagen IV intensity")
    ctrl = comparisons.loc[comparisons.dose_a_um == 0].copy()
    is_mw = method == "mann-whitney"
    method_name = {"mann-whitney": "Mann-Whitney U with Holm correction", "tukey": "Tukey HSD (Tukey-Kramer for unequal group sizes)", "games-howell": "Games-Howell"}[method]
    endpoint = (f"Each point is one gland's log2(day-{int(data.final_day.iloc[0])} bud count / day-{int(data.baseline_day.iloc[0])} bud count), paired by gland identity. "
                "Tests use the displayed log2 ratios, not the untransformed ratios tested in the notebook."
                if kind == "bud" else
                "Each point is the mean raw-image collagen IV AUC across retained line scans from one gland, "
                "divided by the mean of the control gland means. Vehicle mean is 1. "
                "Statistics use these gland means, not individual line scans.")
    family_note = (f"{method_name}: two-sided comparisons of each nonzero dose with DMSO (0). "
                   f"Holm correction covers these {len(comparisons)} comparisons within this plot "
                   "(family-wise alpha 0.05). All tested comparisons are displayed; comparisons between treated doses are not included."
                   if is_mw else
                   f"{method_name}: two-sided all-pairs comparisons across {len(summary)} doses "
                   f"({len(comparisons)} comparisons; family-wise alpha 0.05). The plot shows the {len(summary)-1} "
                   "vehicle comparisons. No second multiplicity correction is applied.")
    notes = [title, endpoint, "Error bar is SD.", family_note,
             "Asterisks: * adjusted p < .05; ** < .01; *** < .001; **** < .0001; ns >= .05. "
             + ("Means and SD are descriptive summaries; the test compares ranks." if is_mw else "Confidence intervals are simultaneous 95% pairwise intervals (approximate for Games-Howell)."),
             "Different glands are treated as independent biological observations. Treatment groups are not paired. "
             "Embryo/litter identifiers are not available for modeling further clustering."]
    if kind == "bud":
        notes.append(f"The {data.cohort.nunique()} source cohorts are pooled as in the notebook. "
                     "Cohort is retained in the exported gland table but is not a covariate in the tests. "
                     "Only the selected baseline/final response is analyzed. For the supplied data this is day 1/day 0; future-washout day-2 counts are excluded.")
    else:
        notes.append("Control normalization uses equal gland weights. The notebook used a pooled control-profile "
                     "mean; the new denominator changes the scale slightly but cannot change these scale-invariant tests.")
    if is_mw:
        notes.extend([diagnostics["method_selection"],
                      "Mann-Whitney does not require normality. It compares independent groups using ranks; "
                      "a pure location/median interpretation additionally requires comparable distribution shapes. "
                      "Unequal variances are not interpreted as evidence specifically about mean differences.",
                      "U is reported with dose B as the first sample. " + diagnostics["rank_effect"]
                      + " Mean differences are descriptive only. No parametric confidence intervals are attached to the rank-test results."])
    else:
        notes.append("Tukey assumes independent approximately normal within-group errors with a common variance. "
                     "Games-Howell allows unequal variances; both require independent glands and approximately normal group means/errors.")
        bf_p = diagnostics["brown_forsythe_p"]
        notes.append(f"Variance diagnostic (median-centered Levene / Brown-Forsythe): "
                     f"F = {diagnostics['brown_forsythe_statistic']:.5g}, p = {bf_p:.5g}; "
                     f"largest/smallest sample variance = {diagnostics['largest_to_smallest_variance']:.3g}. "
                     + diagnostics["small_sample_note"])
        if method == "tukey" and bf_p < .05:
            notes.append("CHECK: The variance diagnostic flags heterogeneity. Interpret the common-variance assumption "
                         "with care; use --test games-howell to inspect an unequal-variance analysis.")
    if kind == "bud":
        excluded = audit.get("excluded_measurements", [])
        placeholders = {x["gland_id"] for x in excluded if x.get("placeholder")}
        notes.append(f"Explicit exclusions: {len(placeholders)} filename-marked placeholder glands; "
                     "no real gland is excluded by a count threshold. The input audit records all excluded observations.")
    else:
        norm = audit["normalization"]
        notes.append(f"{audit['n_scans']} line scans from {len(data)} glands; "
                     f"{audit['minimum_scans_per_gland']}–{audit['maximum_scans_per_gland']} scans per gland. "
                     f"Control denominator = {norm['control_gland_mean_auc']:.10g} AUC units; "
                     f"notebook denominator = {norm['notebook_pooled_control_scan_mean_auc']:.10g}. "
                     f"The normalized scale changes by {100*(norm['new_value_divided_by_notebook_value']-1):.4f}%.")
    notes.append("Input audit and explicit exclusions are recorded in the run JSON. Source notebooks and measurement files are unchanged.")
    columns = (["dose_b_um", "n_a", "n_b", "U_b", "rank_biserial_b_minus_a", "p_raw", "p_adjusted", "p_method", "annotation"] if is_mw else
               ["dose_b_um", "n_a", "n_b", "difference_b_minus_a", "ci95_low", "ci95_high", "p_adjusted", "annotation"])
    display = ctrl[columns]
    text_report = "\n\n".join(notes) + "\n\nGROUP SUMMARY\n" + summary.to_string(index=False) + "\n\nVERSUS VEHICLE\n" + display.to_string(index=False, float_format=lambda x: f"{x:.8g}") + "\n"
    if not alternative.empty:
        alt_ctrl = alternative.loc[alternative.dose_a_um == 0]
        sensitivity_note = ("Sensitivity analysis only (not used for plot asterisks): " + str(alternative.method.iloc[0])
                            + ", using the same gland endpoints and all-pairs comparison family.")
        text_report += "\n" + sensitivity_note + "\n" + alt_ctrl.to_string(index=False, float_format=lambda x: f"{x:.8g}") + "\n"
    if not cohort_sensitivity.empty:
        text_report += "\nCOHORT SENSITIVITY (not used for plot asterisks)\n" + json.dumps(cohort_metadata, indent=2) + "\n" + cohort_sensitivity.to_string(index=False, float_format=lambda x: f"{x:.8g}") + "\n"
    (folder/(prefix+"-statistics.txt")).write_text(text_report)
    css = "body{font:16px/1.55 system-ui,sans-serif;max-width:1000px;margin:36px auto;padding:0 24px;color:#222}table{border-collapse:collapse;font-size:14px}td,th{padding:7px 12px;border-bottom:1px solid #ddd;text-align:right}th{background:#f4f4f4}img{width:410px;max-width:100%}p{max-width:90ch}details{margin-top:24px}"
    image_tag = f'<img src="{html.escape(prefix)}.png" alt="Mean and SD with individual gland points">' if "png" in formats else ""
    html_report = '<!doctype html><html lang="en"><meta charset="utf-8"><title>'+html.escape(title)+'</title><style>'+css+'</style><body><h1>'+html.escape(title)+'</h1>'+image_tag
    html_report += "".join("<p>"+html.escape(x)+"</p>" for x in notes[1:])
    html_report += "<h2>Group summary</h2>"+summary.to_html(index=False, float_format=lambda x: f"{x:.5g}")
    html_report += "<h2>Comparisons against vehicle</h2>"+display.to_html(index=False, float_format=lambda x: f"{x:.6g}")
    detail_label = "Full results: doses versus DMSO (0)" if is_mw else "All pairwise comparisons"
    html_report += "<details><summary>"+detail_label+"</summary>"+comparisons.to_html(index=False, float_format=lambda x: f"{x:.6g}")+"</details>"
    if not alternative.empty:
        html_report += "<h2>Sensitivity analysis</h2><p>"+html.escape(sensitivity_note)+"</p>"+alt_ctrl.to_html(index=False, float_format=lambda x: f"{x:.6g}")
    if not cohort_sensitivity.empty:
        html_report += "<h2>Cohort sensitivity</h2><p>Not used for plot asterisks. Dose and cohort fixed effects, HC3 covariance, t-based inference; Holm correction over the six vehicle comparisons. Contrast confidence intervals here are unadjusted 95%.</p>"+cohort_sensitivity.to_html(index=False, float_format=lambda x: f"{x:.6g}")
    if is_mw:
        html_report += '<p>Method references: <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html">SciPy Mann-Whitney U</a>; <a href="https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html">Holm correction</a>.</p>'
    else:
        html_report += '<p>Method references: <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.tukey_hsd.html">SciPy Tukey / Games-Howell</a>.</p>'
    html_report += '</body></html>\n'
    (folder/(prefix+"-statistics.html")).write_text(html_report)
    data.to_csv(folder/(prefix+"-glands.csv"), index=False)
    summary.to_csv(folder/(prefix+"-summary.csv"), index=False)
    comparisons.to_csv(folder/(prefix+"-comparisons.csv"), index=False)
    ctrl.to_csv(folder/(prefix+"-versus-control.csv"), index=False)
    if not alternative.empty:
        alternative.to_csv(folder/(prefix+"-sensitivity-all-pairs.csv"), index=False)
    if not cohort_sensitivity.empty:
        cohort_sensitivity.to_csv(folder/(prefix+"-cohort-sensitivity.csv"), index=False)


def audit_source_paths(audit):
    """Collect source paths for the publication collision guard."""
    paths = set()
    if isinstance(audit, dict):
        for key, value in audit.items():
            if key in {"path", "source_file", "record_file", "count_file", "config_path"} and isinstance(value, str):
                paths.add(Path(value).resolve())
            else:
                paths.update(audit_source_paths(value))
    elif isinstance(audit, list):
        for value in audit:
            paths.update(audit_source_paths(value))
    return paths


def cli(kind, loader, default_config):
    parser = argparse.ArgumentParser(description="Plot BB-94 titration with gland-level means, SD and adjusted significance; export readable statistics.")
    parser.add_argument("--config", type=Path, default=default_config, help="JSON input manifest; relative input paths resolve beside the manifest.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="Fig4-bud-count-titration" if kind == "bud" else "Fig4-collagen-IV-intensity-titration")
    parser.add_argument("--test", choices=["mann-whitney", "tukey", "games-howell"], default="mann-whitney", help="Default: Mann-Whitney U versus DMSO (0), Holm across those comparisons. Legacy parametric modes retain all-pairs families.")
    parser.add_argument("--formats", nargs="+", choices=["pdf", "svg", "png"], default=["pdf", "svg", "png"])
    parser.add_argument("--font", default="Arial")
    parser.add_argument("--width", type=float, default=1.9, help="Figure width in inches.")
    parser.add_argument("--height", type=float, default=1.75, help="Figure height in inches.")
    parser.add_argument("--overwrite", action="store_true", help="Replace this prefix's generated outputs only.")
    args = parser.parse_args()
    try:
        if not args.prefix or Path(args.prefix).name != args.prefix or args.prefix in {".", ".."}:
            raise ValueError("--prefix must be a plain filename stem, without a directory.")
        if not np.isfinite([args.width, args.height]).all() or min(args.width, args.height) < 1:
            raise ValueError("Figure dimensions must be at least 1 inch.")
        output = args.output_dir.resolve()
        suffixes = ["."+x for x in args.formats] + ["-statistics.txt", "-statistics.html", "-glands.csv", "-summary.csv", "-comparisons.csv", "-all-pairs.csv", "-versus-control.csv", "-sensitivity-all-pairs.csv", "-run.json"]
        if kind == "bud":
            suffixes.append("-cohort-sensitivity.csv")
        targets = [output/(args.prefix+s) for s in suffixes]
        existing = [p for p in targets if p.exists()]
        if existing and not args.overwrite:
            raise ValueError(f"Outputs already exist (e.g. {existing[0]}); choose another folder/prefix or pass --overwrite.")
        data, audit = loader(args.config.resolve())
        sources = audit_source_paths(audit) | {args.config.resolve()}
        if any(p.resolve() in sources for p in targets):
            raise ValueError("An output path collides with an input file. Choose another output folder or prefix.")
        summary, comparisons, diagnostics = analyze(data, args.test)
        alternative_method, alternative = None, pd.DataFrame()
        if args.test != "mann-whitney":
            alternative_method = "tukey" if args.test == "games-howell" else "games-howell"
            _, alternative, _ = analyze(data, alternative_method)
        cohort_sensitivity, cohort_metadata = pd.DataFrame(), {}
        if kind == "bud" and args.test != "mann-whitney":
            from titration_sensitivity import bud_cohort_sensitivity
            cohort_sensitivity, cohort_metadata = bud_cohort_sensitivity(data)
        output.mkdir(parents=True, exist_ok=True)
        # Publish only after every figure and report has been generated successfully.
        with tempfile.TemporaryDirectory(prefix=".titration-", dir=output) as work:
            folder = Path(work)
            plot(data, summary, comparisons, kind, folder/args.prefix, args.formats, args.font, args.width, args.height)
            reports(folder, args.prefix, kind, data, summary, comparisons, diagnostics, audit, args.test, args.formats, alternative, cohort_sensitivity, cohort_metadata)
            run = {"command": sys.argv, "python": platform.python_version(),
                   "packages": {name: importlib.metadata.version(name) for name in ["numpy", "pandas", "scipy", "matplotlib", "seaborn", "statsmodels"]},
                   "config": str(args.config.resolve()), "config_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
                   "analysis": args.test, "alpha": .05, "pairwise_family_size": len(comparisons),
                   "comparison_family": "doses versus DMSO (0)" if args.test == "mann-whitney" else "all pairwise dose comparisons",
                   "endpoint": audit.get("endpoint", "normalized per-gland mean AUC"),
                   "error_bar": "SD", "glands_total": len(data), "diagnostics": diagnostics, "input_audit": audit,
                   "sensitivity_method": alternative_method, "cohort_sensitivity": cohort_metadata,
                   "code_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob("*titration*.py")}}
            (folder/(args.prefix+"-run.json")).write_text(json.dumps(run, indent=2, allow_nan=False, default=str)+"\n")
            for p in folder.iterdir():
                p.replace(output/p.name)
            old_comparison_name = output/(args.prefix+"-all-pairs.csv")
            if old_comparison_name.exists():
                old_comparison_name.unlink()
            if args.test == "mann-whitney":
                # These prior generated supplements do not describe the current run.
                stale_suffixes = ["-sensitivity-all-pairs.csv"]
                if kind == "bud":
                    stale_suffixes.append("-cohort-sensitivity.csv")
                for suffix in stale_suffixes:
                    stale = output/(args.prefix+suffix)
                    if stale.exists():
                        stale.unlink()
        family = "dose-versus-DMSO" if args.test == "mann-whitney" else "all-pairs"
        print(f"Saved {len(data)} glands; {args.test}, {len(comparisons)} {family} comparisons; error bar is SD.")
        print(comparisons.loc[comparisons.dose_a_um == 0, ["dose_b_um", "p_adjusted", "annotation"]].to_string(index=False))
        if "brown_forsythe_p" in diagnostics:
            print(f"Brown-Forsythe variance diagnostic: p={diagnostics['brown_forsythe_p']:.6g}")
        print(f"Report: {output/(args.prefix+'-statistics.html')}")
    except (ValueError, OSError, KeyError) as exc:
        parser.error(str(exc))
