#!/usr/bin/env python3
"""Check bundled sources, accepted experimental analyses and model exports."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import plot_style as style
from scipy.stats import spearmanr, t as student_t
from statsmodels.stats.multitest import multipletests
from dunnett_statistics import _adjusted_tail
ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-models", action="store_true")
    args = parser.parse_args()
    sources = json.loads((ROOT/"source-data-manifest.json").read_text())
    for source in sources:
        path = ROOT/source["path"]
        assert path.is_file(), f"Missing source: {source['path']}"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == source["sha256"], source["path"]
    reference = json.loads((ROOT/"validation/expected-results.json").read_text())
    for spec in reference["tables"]:
        expected = pd.DataFrame(spec["records"])[spec["columns"]].sort_values(spec["sort_by"]).reset_index(drop=True)
        observed = pd.read_csv(ROOT/spec["path"], dtype=spec.get("dtype", {}))[spec["columns"]].sort_values(spec["sort_by"]).reset_index(drop=True)
        pd.testing.assert_frame_equal(observed, expected, check_dtype=False, check_exact=False,
                                      rtol=1e-10, atol=1e-12, obj=spec["path"])
    # Frozen fixtures verify endpoints separately from the new inference.
    dunnett_panels = [("Fig4F", "BB-94-titration/results/bud-count/Fig4F-bud-count-titration", "dose_um", "value"),
                     ("Fig4G", "BB-94-titration/results/collagen-IV/Fig4G-collagen-IV-intensity-titration", "dose_um", "value")]
    if not args.skip_models:
        dunnett_panels.append(("Fig4E", "model-morphometry/results/Fig4E-model-bud-count", "z", "N_buds"))
    for panel, stem, x, y in dunnett_panels:
        suffix = "-replicates.csv" if panel == "Fig4E" else "-glands.csv"
        data = pd.read_csv(ROOT/(stem+suffix))
        groups = [g[y].to_numpy() for _, g in data.groupby(x, sort=True)]
        ns = tuple(len(g) for g in groups)
        df = sum(ns)-len(ns)
        variance = sum(np.sum((g-g.mean())**2) for g in groups)/df
        t = np.array([(g.mean()-groups[0].mean())/np.sqrt(variance*(1/len(g)+1/ns[0])) for g in groups[1:]])
        raw = 2*student_t.sf(np.abs(t), df)
        result = pd.read_csv(ROOT/(stem+"-comparisons.csv"))
        np.testing.assert_allclose(result["t_comparison" if panel == "Fig4E" else "t_b"], t, rtol=1e-12)
        np.testing.assert_allclose(result.p_raw, raw, rtol=1e-10, atol=1e-25)
        assert len(result) == len(groups)-1 and result.df.eq(df).all()
        adjusted = result.p_adjusted.to_numpy()
        assert np.all(adjusted >= raw*(1-1e-7)) and np.all(adjusted <= np.minimum(1, len(result)*raw)*(1+1e-7))
        # Independent quadrature resolution checks the tiny adjusted tails.
        coarse = [_adjusted_tail(float(abs(v)), ns[1:], ns[0], df, nodes=128) for v in t]
        np.testing.assert_allclose(adjusted, coarse, rtol=1e-7, atol=1e-22)
    washout = pd.read_csv(ROOT/"BB-94-washout/results/combined/Fig5C-washout-log2-ratios-data.csv")
    expected_welch = []
    for stage in ["first", "second"]:
        a = washout.loc[(washout.stage == stage) & (washout.group == "ctrl"), "log2_ratio"].to_numpy()
        b = washout.loc[(washout.stage == stage) & (washout.group == "BB94"), "log2_ratio"].to_numpy()
        va, vb = a.var(ddof=1)/len(a), b.var(ddof=1)/len(b)
        t = (b.mean()-a.mean())/np.sqrt(va+vb)
        df = (va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1))
        expected_welch.append([t, df, 2*student_t.sf(abs(t), df)])
    expected_welch = np.array(expected_welch)
    for variant in ["uniform-stars", "symbols-stars-no-legend"]:
        path = ROOT/"BB-94-washout/results/combined"/("Fig5C-washout-log2-ratios-"+variant+"-statistics.csv")
        actual = pd.read_csv(path)
        np.testing.assert_allclose(actual[["t_BB94_minus_control", "df", "p_raw"]], expected_welch, rtol=1e-12)
        np.testing.assert_allclose(actual.p_Holm, multipletests(expected_welch[:, 2], method="holm")[1], rtol=1e-12)
    native_plots = []
    for palette in style.palettes()["titration"]:
        suffix = "-"+palette
        native_plots.append(("BB-94-titration/results/combined/Fig4FG-BB94-titration"+suffix,
                             "fig4-titration-combined", palette))
    for name in ["uniform", "symbols"]:
        stem = ("Fig5C-washout-log2-ratios-uniform-stars" if name == "uniform"
                else "Fig5C-washout-log2-ratios-symbols-stars-no-legend")
        native_plots.append(("BB-94-washout/results/combined/"+stem, "fig5-washout-"+name, None))
    for stem, key, palette in native_plots:
        for extension in ["pdf", "svg", "png"]:
            path = ROOT/(stem+"."+extension)
            assert path.is_file() and path.stat().st_size > 0, path
        path = ROOT/(stem+".svg")
        svg = ET.parse(path).getroot()
        colors = (style.titration_colors(palette, style.palettes()["titration_dose_order_um"])
                  if palette else [style.palettes()["washout"][group]
                                   for group in ["control", "BB94"]])
        for color in colors:
            assert color.lower() in path.read_text().lower(), (stem, color)
        size = [float(svg.get(k).removesuffix("pt")) for k in ["width", "height"]]
        np.testing.assert_allclose(size, np.array(style.dimensions(key))*72, atol=1e-6)
        if palette:
            annotations = {n.get("id"): n.get("aria-label") for n in svg.iter()
                           if n.get("id") in {"compact-x-prefix", "compact-x-unit"}}
            assert annotations == {"compact-x-prefix": "[BB-94]", "compact-x-unit": "(µM)"}
            assert svg.get("overflow") == "visible"
            groups = [n for n in svg.iter() if n.get("id", "").startswith("PathCollection_")]
            counts = [sum(n.tag.endswith("}use") for n in group.iter()) for group in groups]
            assert len(counts) == 14 and sum(counts[:7]) == 110 and sum(counts[7:]) == 49
    combined = json.loads((ROOT/"BB-94-titration/results/combined/Fig4FG-BB94-titration-run.json").read_text())
    for path, digest in combined["input_table_sha256"].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    subprocess.run([sys.executable, str(ROOT/"mesh-morphometry/verify_mesh_plots.py")], cwd=ROOT, check=True)
    mesh = json.loads((ROOT/"mesh-morphometry/validation/last-verification.json").read_text())
    model = {"plots_verified": 0}
    if not args.skip_models:
        subprocess.run([sys.executable, str(ROOT/"model-morphometry/verify_model_plots.py")], cwd=ROOT, check=True)
        model = json.loads((ROOT/"model-morphometry/validation/last-verification.json").read_text())
    report = {"status": "passed", "source_files_verified": len(sources),
              "dunnett_panels_verified": [p[0] for p in dunnett_panels],
              "welch_fig5C_comparisons_verified": 2,
              "numerical_tables_verified": len(reference["tables"]),
              "plots_verified": len(native_plots)+mesh["plots_verified"]+model["plots_verified"],
              "titration_color_styles_verified": list(style.palettes()["titration"]),
              "mesh_curves_verified": mesh["individual_and_mean_curves_verified"],
              "mesh_sd_bands_verified": mesh["sd_bands_verified"],
              "mesh_recolored_curves_verified": mesh["recolored_curves_verified"],
              "mesh_recolored_sd_bands_verified": mesh["recolored_sd_bands_verified"],
              "model_verification": model, "combined_titration_gland_points": [110, 49],
              "compact_markers": {"shape": "hollow circle", "diameter_points": {"4E": style.COMPACT_CIRCLE_DIAMETER,
                                                                                  "4F": style.COMPACT_BUD_RATIO_CIRCLE_DIAMETER,
                                                                                  "4G": style.COMPACT_CIRCLE_DIAMETER},
                                   "edge_width_points": style.COMPACT_CIRCLE_EDGE_WIDTH,
                                   "summary_color": "black", "summary_above_observations": True},
              "relative_tolerance": 1e-10, "absolute_tolerance": 1e-12}
    trend_report = json.loads((ROOT/"validation/Fig4EFG-spearman-trends.json").read_text())
    trend_tables = {
        "Fig4E": ("model-morphometry/results/Fig4E-model-bud-count-replicates.csv", "z", "N_buds"),
        "Fig4F": ("BB-94-titration/results/bud-count/Fig4F-bud-count-titration-glands.csv", "dose_um", "value"),
        "Fig4G": ("BB-94-titration/results/collagen-IV/Fig4G-collagen-IV-intensity-titration-glands.csv", "dose_um", "value"),
    }
    panel_order = ["Fig4E", "Fig4F", "Fig4G"]
    adjusted = multipletests([trend_report["panels"].get(p, {}).get("p_raw", 1) for p in panel_order], method="holm")[1]
    for i, panel in enumerate(panel_order):
        if args.skip_models and panel == "Fig4E":
            continue
        trend = trend_report["panels"][panel]
        path, x, y = trend_tables[panel]
        data = pd.read_csv(ROOT/path)
        assert trend["n"] == len(data)
        np.testing.assert_allclose(trend["rho"], spearmanr(data[x], data[y]).statistic, atol=1e-12)
        assert trend["permutations"] == 999999 and trend["permutation_mode"] == "Monte Carlo"
        assert trend["seed"] == 20261018
        np.testing.assert_allclose(trend["p_raw"], (trend["extreme_permutations"]+1)/(trend["permutations"]+1))
        np.testing.assert_allclose(trend["p_adjusted"], adjusted[i])
    report["figure4_spearman_trends_verified"] = [p for p in panel_order if p in trend_report["panels"]]
    vector_check = [sys.executable, str(ROOT/"verify_vector_exports.py")]
    if args.skip_models:
        vector_check.append("--skip-models")
    subprocess.run(vector_check, cwd=ROOT, check=True)
    vector_report = json.loads((ROOT/"validation/vector-outline-checks.json").read_text())
    report["outlined_pdf_svg_pairs_verified"] = vector_report["plot_sets_verified"]
    report["export_font_resources"] = vector_report["font_resources"]
    report["export_raster_images"] = vector_report["raster_images"]
    (ROOT/"validation/last-verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(f"Verified {len(sources)} source files, {len(reference['tables'])} experimental tables and {report['plots_verified']} plot exports.")


if __name__ == "__main__":
    main()
