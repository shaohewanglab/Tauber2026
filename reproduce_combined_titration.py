#!/usr/bin/env python3
"""Stack experimental 4F/G with a shared concentration axis and no x title."""
import hashlib
import argparse
import json
from pathlib import Path
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/"BB-94-titration/plotting_scripts"))
from titration_analysis import draw_panel
import plot_style as style
from trend_statistics import permutation_spearman


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT/"BB-94-titration/results/combined")
    parser.add_argument("--font", default="Arial")
    args = parser.parse_args()
    results = args.output_dir.resolve()
    results.mkdir(parents=True, exist_ok=True)
    panels, inputs = [], {}
    for kind, folder, prefix in [("bud", "bud-count", "Fig4F-bud-count-titration"),
                                 ("collagen", "collagen-IV", "Fig4G-collagen-IV-intensity-titration")]:
        tables = []
        for suffix in ["glands", "summary", "comparisons"]:
            path = ROOT/"BB-94-titration/results"/folder/(prefix+"-"+suffix+".csv")
            inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
            tables.append(pd.read_csv(path))
        panels.append((kind, *tables))
    np.testing.assert_array_equal(panels[0][2].dose_um, panels[1][2].dose_um)
    trends = {}
    for kind, data, summary, comparisons in panels:
        panel = "Fig4F" if kind == "bud" else "Fig4G"
        result = permutation_spearman(data.dose_um, data.value)
        result.update(panel=panel, predictor="BB-94 concentration (uM)",
                      outcome="log2(day1/day0 bud-count ratio)" if kind == "bud" else "Normalized per-gland collagen IV intensity",
                      statistical_unit="One gland", input_table_sha256=inputs)
        trends[panel] = result
    (results/"Fig4FG-BB94-titration-trends.json").write_text(json.dumps(trends, indent=2)+"\n")
    for palette in style.palettes()["titration"]:
        with mpl.rc_context(style.rc(font=args.font)):
            fig, axes = plt.subplots(2, 1, sharex=True,
                                     figsize=style.dimensions("fig4-titration-combined"), dpi=300)
            for ax, (kind, data, summary, comparisons) in zip(axes, panels):
                draw_panel(ax, data, summary, comparisons, kind, palette)
            style.stacked_titration_layout(fig, axes)
            fig.canvas.draw()
            for ax, (_, data, summary, _) in zip(axes, panels):
                points = [c for c in ax.collections if hasattr(c, "get_offsets") and len(c.get_offsets()) > 1]
                assert len(points) == 7
                for dose, collection in zip(summary.dose_um, points):
                    np.testing.assert_allclose(collection.get_offsets()[:, 1], data.loc[data.dose_um == dose, "value"],
                                               rtol=1e-12, atol=1e-12)
            suffix = "-"+palette
            stem = results/("Fig4FG-BB94-titration"+suffix)
            style.save_native_figure(fig, stem, dpi=450, allow_overflow=True)
            plt.close(fig)
    record = {"panels": ["4F: log2 bud ratio", "4G: collagen IV intensity"],
              "full_figure_inches": style.dimensions("fig4-titration-combined"),
              "x_axis_shared": True, "x_axis_titles": False,
              "x_tick_labels": "Concentrations shown on lower panel only",
              "x_tick_side_annotations": {"left_lines": ["[BB-94]", "(µM)"], "font_points": 6,
                                          "placement": "Both lines right-aligned to left edge of y tick labels; canvas overflow allowed"},
              "markers": {"shape": "hollow circle", "diameter_points": {"4F": style.COMPACT_BUD_RATIO_CIRCLE_DIAMETER,
                                                                         "4G": style.COMPACT_CIRCLE_DIAMETER},
                          "edge_width_points": style.COMPACT_CIRCLE_EDGE_WIDTH,
                          "summary_color": "black", "summary_above_observations": True},
              "font": args.font, "axis_font_points": 6,
              "export_text": "Vector outlines in PDF and SVG; no fonts required to view or place",
              "data_and_statistics": "Gland endpoints unchanged; two-sided single-step Dunnett tests versus vehicle, separate six-comparison families for 4F and 4G",
              "input_table_sha256": inputs, "gland_points_per_panel": [110, 49],
              "color_styles": list(style.palettes()["titration"])}
    (results/"Fig4FG-BB94-titration-run.json").write_text(json.dumps(record, indent=2)+"\n")
    print("Saved the Tol Bright shared-axis 4F/G plot as PDF/SVG/PNG; verified all 159 gland coordinates.")


if __name__ == "__main__":
    main()
