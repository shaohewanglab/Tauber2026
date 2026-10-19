#!/usr/bin/env python3
"""Verify model sources, frozen original curves, bud counts and exports."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
from collections import Counter
from matplotlib.colors import to_hex
import numpy as np
import pandas as pd
import reproduce_model_plots as plots

ROOT = Path(__file__).resolve().parent


def horizontal_frame_bounds(path):
    """Read the two long horizontal spines from the outlined TeX SVG."""
    bounds = []
    for node in ET.parse(path).getroot().iter():
        if not node.tag.endswith("}path") or node.get("id") or node.get("fill") != "none":
            continue
        match = re.fullmatch(r"M([-\d.]+) [-\d.]+H([-\d.]+)", node.get("d", ""))
        if match:
            bounds.append(tuple(map(float, match.groups())))
    pairs = [pair for pair, count in Counter(bounds).items() if count == 2]
    assert pairs, f"No paired horizontal frame spines: {path}"
    return max(pairs, key=lambda pair: pair[1]-pair[0])


def main():
    sources = json.loads((ROOT/"source-data-manifest.json").read_text())
    for source in sources:
        assert hashlib.sha256((ROOT/source["path"]).read_bytes()).hexdigest() == source["sha256"], source["path"]
    references = json.loads((ROOT/"validation/original-curves.json").read_text())
    curves = pd.read_csv(ROOT/"results/acircularity-curves.csv")
    for panel, source in references.items():
        assert hashlib.sha256((ROOT/"original_scripts"/source["source_script"]).read_bytes()).hexdigest() == source["script_sha256"]
        groups = list(curves.loc[curves.panel == panel].groupby(["condition", "z"], sort=False))
        assert len(groups) == len(source["lines"])
        for (_, actual), original in zip(groups, source["lines"], strict=True):
            np.testing.assert_allclose(actual.tau_relative, original["x"], rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(actual.delta_chi_percent, original["y"], rtol=1e-12, atol=1e-12)
    buds, summary = plots.model_bud_data()
    pd.testing.assert_frame_equal(pd.read_csv(ROOT/"results/Fig4E-model-bud-count-replicates.csv"),
                                  buds, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    pd.testing.assert_frame_equal(pd.read_csv(ROOT/"results/Fig4E-model-bud-count-summary.csv"),
                                  summary, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    comparisons = plots.analyze_bud_counts(buds, plots.Z_VALUES)
    pd.testing.assert_frame_equal(pd.read_csv(ROOT/"results/Fig4E-model-bud-count-comparisons.csv",
                                             keep_default_na=False, na_values=[""]),
                                  comparisons, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    statistics = json.loads((ROOT/"results/Fig4E-model-bud-count-statistics.json").read_text())
    assert statistics["source_sha256"] == hashlib.sha256((ROOT/"ModelData/bud_counts_summary.csv").read_bytes()).hexdigest()
    assert statistics["baseline_z"] == min(plots.Z_VALUES) and statistics["family_size"] == 5
    assert statistics["comparisons_tested"] == int((comparisons.status == "tested").sum())
    assert statistics["unavailable_comparisons"] == int((comparisons.status != "tested").sum())
    for row, actual in zip(comparisons.itertuples(), statistics["results"], strict=True):
        assert row.annotation == actual["annotation"]
        if row.status != "tested":
            assert actual["p_raw"] is None and actual["p_adjusted"] is None
        else:
            np.testing.assert_allclose([actual["p_raw"], actual["p_adjusted"]],
                                       [row.p_raw, row.p_adjusted], rtol=1e-12, atol=1e-12)
    exports = [(plots.STEMS["4A"]+"-tol-bright", "fig4-model-2d"),
               (plots.STEMS["5D"], "fig5-future-d"), (plots.STEMS["5E"], "fig5-future-e")]
    for palette in plots.style.palettes()["titration"]:
        exports.append(("Fig4E-model-bud-count-"+palette, "fig4-model-buds"))
    for stem, key in exports:
        path = ROOT/"results"/stem
        for extension in ["pdf", "svg", "png"]:
            assert path.with_suffix("."+extension).is_file(), path
        svg = ET.parse(path.with_suffix(".svg")).getroot()
        actual = [float(svg.get(dim).removesuffix("pt"))/72 for dim in ["width", "height"]]
        np.testing.assert_allclose(actual, plots.style.dimensions(key), rtol=0, atol=1e-4)
        if "tol-bright" in stem:
            text = re.sub(r"#[0-9a-fA-F]{3,8}\b", lambda m: to_hex(m.group()).upper(),
                          path.with_suffix(".svg").read_text())
            colors = plots.style.thickness_colors("tol-bright", plots.Z_VALUES)
            positions = [text.index(color) for color in colors]
            assert positions == sorted(positions), stem
    alignment = {}
    for suffix in ["-tol-bright"]:
        bounds = horizontal_frame_bounds(ROOT/"results"/(plots.STEMS["4A"]+suffix+".svg"))
        width_inches = (bounds[1]-bounds[0])/72.27  # dvisvgm's PGF coordinates use TeX points.
        # TeX quantizes physical coordinates; compare the configured nominal
        # width within 0.00001 inch, then require exact exported equality to 4B.
        np.testing.assert_allclose(width_inches, plots.style.frame_width_inches("fig4-mesh"),
                                   rtol=0, atol=1e-5)
        reference = ROOT.parent/"mesh-morphometry/results"/("Fig4B-asphericitychange_evolution_BB94_cleanmesh_plus"+suffix+".svg")
        if reference.is_file():
            np.testing.assert_array_equal(bounds, horizontal_frame_bounds(reference))
        alignment[suffix.removeprefix("-")] = {"frame_width_inches": width_inches,
                                            "horizontal_bounds_tex_points": bounds,
                                            "matches_exported_4B": reference.is_file()}
    report = {"status": "passed", "source_files_verified": len(sources),
              "original_model_curves_verified": 11, "plots_verified": len(exports),
              "bud_counts_available": len(buds), "bud_counts_target": 3*len(summary),
              "bud_counts_missing": 3*len(summary)-len(buds),
              "missing_sd_at_z": summary.loc[summary.V_count == 1, "z"].tolist(),
              "internal_data_note": "INTERNAL_NOTES.md",
              "fig4A_frame_alignment": alignment,
              "fig4E_statistics": {"baseline_z": statistics["baseline_z"], "family_size": 5,
                                     "comparisons_tested": statistics["comparisons_tested"],
                                     "unavailable_comparisons": statistics["unavailable_comparisons"]},
              "relative_tolerance": 1e-12, "absolute_tolerance": 1e-12}
    (ROOT/"validation/last-verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(f"Verified 13 model source/script files, 11 original curves, {len(buds)} available bud counts and {len(exports)} model plots.")


if __name__ == "__main__":
    main()
