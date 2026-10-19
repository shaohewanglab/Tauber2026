#!/usr/bin/env python3
"""Reproduce model acircularity and model bud-count plots."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy
import sys
from unittest.mock import patch

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_style as style
from model_bud_stats import analyze_bud_counts, MIN_REPLICATES
from trend_statistics import permutation_spearman
from dunnett_statistics import ASSUMPTIONS

ROOT = Path(__file__).resolve().parent
Z_VALUES = [.0005, .001, .003, .005, .008, .01]
SCRIPTS = {"4A": "create_acircularity_plot_exp.py",
           "5D": "create_acircularity_plot_washout.py",
           "5E": "create_acircularity_plot_washoutcontrol.py"}
STEMS = {"4A": "Fig4A-acircularity-varying-thickness",
         "5D": "Fig5D-acircularity-full-washout",
         "5E": "Fig5E-acircularity-post-washout"}


def source_path(z, washout=False):
    folder = "washout_paper" if washout else "thickness_paper"
    dt = "0.01" if z == .0005 else "0.1"
    suffix = "_washout" if washout else ""
    return ROOT / "ModelData" / folder / (
        f"shape_array_r_5_z_{z}_kP_0.01_kn_5e-06_tg_1e+06_tgi_2e+06_"
        f"M_1e-07_dt_{dt}_expV0.1{suffix}.npz")


def measures(z, washout=False):
    with np.load(source_path(z, washout), allow_pickle=False) as archive:
        shapes = archive["data"]
    if shapes.ndim != 3 or shapes.shape[1] != 2 or not np.isfinite(shapes).all():
        raise ValueError("Expected finite (node, xy, saved-frame) model coordinates")
    area, perimeter = np.zeros(shapes.shape[2]), np.zeros(shapes.shape[2])
    for frame in range(shapes.shape[2]):
        contour = shapes[:, :, frame]
        # Preserve Justin's convention: close the area with a chord, but sum
        # only consecutive-node edges for the perimeter of this open contour.
        x = np.append(contour[:, 0], contour[0, 0])
        y = np.append(contour[:, 1], contour[0, 1])
        area[frame] = .5 * np.abs(np.sum(x[:-1]*y[1:] - y[:-1]*x[1:]))
        perimeter[frame] = np.sum(np.sqrt(np.sum((contour[1:]-contour[:-1])**2, axis=1)))
    tau = np.arange(shapes.shape[2]) * .01
    chi = 1 / ((2*np.sqrt(np.pi*area))/perimeter) - 1
    return tau, area, perimeter, chi


def capture_original_curves():
    """Execute the supplied scripts, redirecting only loads and output calls."""
    original_load = np.load

    def redirected_load(path, *args, **kwargs):
        return original_load(ROOT / "ModelData" / str(path).removeprefix("output/"), *args, **kwargs)

    sys.path.insert(0, str(ROOT / "original_scripts"))
    references = {}
    for panel, script in SCRIPTS.items():
        with mpl.rc_context(), patch.object(np, "load", redirected_load), \
                patch.object(plt, "savefig"), patch.object(plt, "show"):
            namespace = runpy.run_path(str(ROOT / "original_scripts" / script))
            ax = namespace["ax"]
            references[panel] = {
                "source_script": script,
                "script_sha256": hashlib.sha256((ROOT / "original_scripts" / script).read_bytes()).hexdigest(),
                "lines": [{"x": line.get_xdata().tolist(), "y": line.get_ydata().tolist()}
                          for line in ax.lines],
                "bands": [collection.get_paths()[0].vertices.tolist() for collection in ax.collections],
            }
            plt.close(namespace["fig"])
    return references


def add_curve(ax, panel, z, washout, color, label, records, start=0, zorder=2):
    tau, area, perimeter, chi = measures(z, washout)
    x = tau[start:] - tau[start]
    y = 100 * (chi[start:] - chi[start])
    ax.plot(x, y, color=color, label=label, lw=style.MEAN_WIDTH, zorder=zorder)
    records.extend({"panel": panel, "condition": "washout" if washout else "reference",
                    "z": z, "saved_frame": i, "tau_relative": t, "delta_chi_percent": value}
                   for i, t, value in zip(range(start, len(tau)), x, y))


def build_acircularity(panel, palette="tol-bright"):
    mpl.rcParams.update(style.rc(arev_math=True))
    mpl.rcParams.update({"xtick.major.pad": 1, "ytick.major.pad": 1,
                         "axes.labelpad": style.MESH_AXIS_LABEL_PAD})
    key = {"4A": "fig4-model-2d", "5D": "fig5-future-d", "5E": "fig5-future-e"}[panel]
    fig, ax = plt.subplots(figsize=style.dimensions(key))
    records = []
    if panel == "4A":
        colors = style.thickness_colors(palette, Z_VALUES)
        for z, color in zip(Z_VALUES, colors):
            add_curve(ax, panel, z, False, color, f"{z:g}", records)
        ax.legend(title=r"$z$", ncols=2, loc="upper left", frameon=False,
                  handlelength=.6, handletextpad=.35, columnspacing=.7,
                  labelspacing=.25, borderaxespad=0)
    elif panel == "5D":
        add_curve(ax, panel, .015, True, "firebrick", "Washout", records, zorder=3)
        for z, gray in [(.001, .44), (.015, .70)]:
            add_curve(ax, panel, z, False, mpl.colormaps["Greys"](gray), rf"$z$ = {z:g}", records)
        ax.legend(loc="upper left", frameon=False, handlelength=.6, handletextpad=.4)
    else:
        add_curve(ax, panel, .015, True, "firebrick", "Post washout", records, start=75, zorder=3)
        add_curve(ax, panel, .001, False, mpl.colormaps["Greys"](.44), r"$z$ = 0.001", records, start=75)
        ax.legend(loc="upper left", frameon=False, handlelength=.6, handletextpad=.4)
    ax.set_xlabel(r"$\Delta\tau/\tau_g$" if panel == "5E" else r"$\tau/\tau_g$")
    ax.set_ylabel(r"$\Delta\chi_{\text{2D}}$ [%]")
    ax.set_xlim(0, .75 if panel == "5E" else 1.5)
    ax.set_ylim(0, 110 if panel == "4A" else 150)
    for axis, major, minor in [(ax.xaxis, .5, .1), (ax.yaxis, 50, 10)]:
        axis.set_major_locator(MultipleLocator(major))
        axis.set_minor_locator(MultipleLocator(minor))
    ax.tick_params(which="major", length=3, direction="out")
    ax.tick_params(which="minor", length=1.5, direction="out")
    if panel == "4A":
        # The endpoint label otherwise reserves extra right margin, narrowing
        # the frame relative to 4B despite their identical full canvas sizes.
        for value, label in zip(ax.get_xticks(), ax.get_xticklabels()):
            if np.isclose(value, ax.get_xlim()[1]):
                label.set_ha("right")
    style.compact_layout(fig, ax, key, arev_math=True)
    return fig, records


def model_bud_data():
    summary = pd.read_csv(ROOT / "ModelData/bud_counts_summary.csv")
    rows = []
    for row in summary.itertuples():
        counts = [float(v) for v in [row.Nbud1, row.Nbud2, row.Nbud3] if np.isfinite(v)]
        assert len(counts) == row.V_count
        np.testing.assert_allclose(np.mean(counts), row.mean, rtol=1e-12)
        if len(counts) > 1:
            np.testing.assert_allclose(np.std(counts, ddof=1), row.std, rtol=1e-12)
        else:
            assert np.isnan(row.std), "SD must remain unavailable with one replicate"
        rows.extend({"z": row.z, "replicate": i+1, "N_buds": value} for i, value in enumerate(counts))
    return pd.DataFrame(rows), summary


def build_model_buds(palette="tol-bright"):
    import seaborn as sns
    mpl.rcParams.update(style.rc(arev_math=True))
    data, summary = model_bud_data()
    comparisons = analyze_bud_counts(data, Z_VALUES)
    fig, ax = plt.subplots(figsize=style.dimensions("fig4-model-buds"))
    colors = style.thickness_colors(palette, Z_VALUES)
    sns.swarmplot(data=data, x="z", y="N_buds", order=Z_VALUES,
                  color="gray", marker="o", size=style.COMPACT_CIRCLE_DIAMETER,
                  linewidth=style.COMPACT_CIRCLE_EDGE_WIDTH, warn_thresh=0, ax=ax)
    fig.canvas.draw()
    for collection, color in zip(ax.collections, colors):
        style.compact_circle_style(collection, color)
    for i, row in enumerate(summary.itertuples()):
        ax.errorbar(i, row.mean, yerr=row.std if np.isfinite(row.std) else None,
                    fmt="_", color="black", markersize=3, markeredgewidth=style.MEAN_WIDTH,
                    elinewidth=style.ERROR_WIDTH, capsize=1, capthick=style.ERROR_WIDTH,
                    zorder=style.COMPACT_SUMMARY_ZORDER)
    ax.set_xticks(range(6), labels=[f"{z:g}" for z in Z_VALUES], rotation=45, ha="right")
    highest = max(float(data.N_buds.max()), float((summary["mean"]+summary["std"]).max()))
    lowest = min(float(data.N_buds.min()), float((summary["mean"]-summary["std"]).min()))
    span = max(highest-lowest, 1)
    annotation_y = highest+.08*span
    for i, row in enumerate(comparisons.itertuples(), start=1):
        ax.text(i, annotation_y, row.annotation, ha="center", va="bottom",
                fontsize=6*style.TEX_POINT_FACTOR, gid=f"bud-stat-{i}")
    ax.set_ylim(0, max(6.7, annotation_y+.16*span))
    ax.set_yticks([0, 3, 6])
    ax.set_ylabel(r"$N_{\text{buds}}$")
    ax.tick_params(pad=2, length=3)
    ax.spines[["right", "top"]].set_visible(False)
    style.compact_titration_layout(fig, ax, "fig4-model-buds", arev_math=True)
    return fig, data, summary


def write_model_bud_statistics(data, output):
    comparisons = analyze_bud_counts(data, Z_VALUES)
    stem = output/"Fig4E-model-bud-count"
    comparisons.to_csv(stem.with_name(stem.name+"-comparisons.csv"), index=False)
    records = [{key: (None if pd.isna(value) else value.item() if isinstance(value, np.generic) else value)
                for key, value in row.items()} for row in comparisons.to_dict(orient="records")]
    report = {"panel": "Fig4E", "baseline_z": min(Z_VALUES),
              "baseline_simulations": int((data.z == min(Z_VALUES)).sum()),
              "method": "Two-sided single-step Dunnett comparisons of mean bud count; five thicker-z versus thinnest-z comparisons",
              "method_selection": "Deterministic conditional normal/chi-square quadrature; Dunnett-adjusted P values, without additional Holm correction",
              "minimum_replicates_per_group": MIN_REPLICATES, "family_size": len(Z_VALUES)-1,
              "comparisons_tested": int((comparisons.status == "tested").sum()),
              "unavailable_comparisons": int((comparisons.status != "tested").sum()),
              "unavailable_annotation": "n/a; full family not tested until every group has at least two independent simulations",
              "available_simulations": len(data), "intended_simulations": 3*len(Z_VALUES),
              "data_status": "complete" if len(data) == 3*len(Z_VALUES) else "incomplete",
              "simulations_per_thickness": {f"{z:g}": int((data.z == z).sum()) for z in Z_VALUES},
              "count_endpoint_tau_over_tau_g": 1.5,
              "source": "ModelData/bud_counts_summary.csv",
              "source_sha256": hashlib.sha256((ROOT/"ModelData/bud_counts_summary.csv").read_bytes()).hexdigest(),
              "alpha": .05, "annotation_thresholds": {"*": "p < .05", "**": "p < .01", "***": "p < .001", "****": "p < .0001"},
              "interpretation": "Compare mean bud counts of independent simulations. No missing observations are imputed.",
              "assumptions": ASSUMPTIONS,
              "small_sample_note": "Three runs per thickness cannot establish the normal/common-variance approximation. Some groups have zero observed SD; the pooled within-group variance is positive. An unbudded shape counts as 1.",
              "results": records}
    stem.with_name(stem.name+"-statistics.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    text = ["FIGURE 4E: MODEL BUD COUNTS VERSUS THINNEST z", report["method"],
            report["method_selection"], f"Minimum independent simulations per group: {MIN_REPLICATES}.",
            report["unavailable_annotation"], report["assumptions"], report["small_sample_note"], "", comparisons.to_string(index=False)]
    stem.with_name(stem.name+"-statistics.txt").write_text("\n".join(text)+"\n")
    return report


def check_curves(fig, expected):
    assert len(fig.axes[0].lines) == len(expected["lines"])
    for actual, source in zip(fig.axes[0].lines, expected["lines"]):
        np.testing.assert_allclose(actual.get_xdata(), source["x"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(actual.get_ydata(), source["y"], rtol=1e-12, atol=1e-12)
    if expected["source_script"] == SCRIPTS["4A"]:
        assert len(fig.axes[0].collections) == 0, "4A's theoretical shading was removed for presentation"
        return
    assert len(fig.axes[0].collections) == len(expected["bands"])
    for actual, source in zip(fig.axes[0].collections, expected["bands"]):
        np.testing.assert_allclose(actual.get_paths()[0].vertices, source, rtol=1e-12, atol=1e-12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sources = json.loads((ROOT/"source-data-manifest.json").read_text())
    for source in sources:
        assert hashlib.sha256((ROOT/source["path"]).read_bytes()).hexdigest() == source["sha256"]
    # The reference was captured from the unchanged supplied scripts. Preserve
    # it as an independent fixture instead of overwriting it during a rebuild.
    expected = json.loads((ROOT/"validation/original-curves.json").read_text())
    rows = []
    for panel in ["4A", "5D", "5E"]:
        with mpl.rc_context():
            fig, records = build_acircularity(panel)
            check_curves(fig, expected[panel])
            rows.extend(records)
            suffix = "-tol-bright" if panel == "4A" else ""
            style.save_arev_figure(fig, args.output_dir/(STEMS[panel]+suffix), dpi=450)
            plt.close(fig)
    pd.DataFrame(rows).to_csv(args.output_dir/"acircularity-curves.csv", index=False)
    with mpl.rc_context():
        fig, data, summary = build_model_buds()
        style.save_arev_figure(fig, args.output_dir/"Fig4E-model-bud-count-tol-bright", dpi=450)
        plt.close(fig)
    data.to_csv(args.output_dir/"Fig4E-model-bud-count-replicates.csv", index=False)
    summary.to_csv(args.output_dir/"Fig4E-model-bud-count-summary.csv", index=False)
    bud_statistics = write_model_bud_statistics(data, args.output_dir)
    trend = permutation_spearman(data.z, data.N_buds)
    trend.update(panel="Fig4E", predictor="Effective thickness z", outcome="Bud count at tau/tau_g = 1.5",
                 statistical_unit="One independent simulation",
                 source="ModelData/bud_counts_summary.csv", source_sha256=bud_statistics["source_sha256"])
    (args.output_dir/"Fig4E-model-bud-count-trend.json").write_text(json.dumps(trend, indent=2)+"\n")
    report = {"source_hashes_checked": len(sources), "original_curves_checked": 11,
              "tol_thickness_curves_checked": 6, "upper_bound_bands_checked": 0,
              "presentation_changes": {"4A": "No gray shading; y = 0 to 110; two legend columns of three entries, shifted up by 1.8 pt. Tol colors run purple, blue, cyan, green, yellow, red as z increases."},
              "tol_thickness_colors": style.palettes()["thickness"]["tol-bright"],
              "model_bud_replicates": len(data), "model_bud_groups": len(summary),
              "model_bud_target_replicates_per_group": 3,
              "model_bud_missing_replicates": int(3 * len(summary) - len(data)),
              "compact_markers": {"shape": "hollow circle", "diameter_points": style.COMPACT_CIRCLE_DIAMETER,
                                   "edge_width_points": style.COMPACT_CIRCLE_EDGE_WIDTH,
                                   "summary_color": "black", "summary_above_observations": True},
              "missing_sd_at_z": summary.loc[summary.V_count == 1, "z"].tolist(),
              "model_plots_exported": 4,
              "figure4_color_style": "tol-bright",
              "fig4A_frame_width_mm": style.frame_width_inches("fig4-model-2d")*25.4,
              "fig4E_statistics": {"baseline_z": bud_statistics["baseline_z"], "family_size": bud_statistics["family_size"],
                                    "comparisons_tested": bud_statistics["comparisons_tested"],
                                    "unavailable_comparisons": bud_statistics["unavailable_comparisons"]},
              "status": "passed", "original_script_changes_for_reference": "Redirect np.load paths and suppress save/show only; numerical code unchanged."}
    (ROOT/"validation/model-verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print("Verified 11 source curves; exported 4 model plots as PDF/SVG/PNG with Tol Bright for Figure 4.")


if __name__ == "__main__":
    main()
