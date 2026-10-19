#!/usr/bin/env python3
"""Check source hashes, all curves/bands against Justin's scripts, and exports."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import reproduce_mesh_plots as plots
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def check(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--report", type=Path, default=ROOT / "validation/last-verification.json")
    args = parser.parse_args()
    sources = json.loads((ROOT / "source-data-manifest.json").read_text())
    used = {f"source_data/volume_area_{sample}_cp_model0907_v2.txt"
            for cohort in plots.COHORTS.values() for sample in cohort}
    check(used == {s["path"] for s in sources}, "Source manifest differs from plotted cohorts")
    for source in sources:
        path = ROOT / source["path"]
        check(path.is_file(), f"Missing source: {path}")
        check(hashlib.sha256(path.read_bytes()).hexdigest() == source["sha256"],
              f"Source changed: {path}")
    reference = json.loads((ROOT / "validation/expected-curves.json").read_text())
    tol_colors = plots.style.thickness_colors("tol-bright", [.001, .005, .01])
    line_count = band_count = summary_count = 0
    for number, expected in reference.items():
        with mpl.rc_context():
            # Numerical verification needs neither TeX nor an installed Arial.
            # Final glyph outlines are audited separately in validation/.
            fig, rows = plots.build_figure(number, "mathtext", layout=False,
                                          font="DejaVu Sans")
            ax = fig.axes[0]
            for observed, wanted in [
                (fig.get_size_inches(), plots.style.dimensions(f"fig{number}-mesh")),
                (ax.get_xlim(), expected["xlim"]),
                (ax.get_ylim(), plots.PRESENTATION_YLIMS.get(number, expected["ylim"])),
            ]:
                np.testing.assert_allclose(observed, wanted, rtol=1e-10, atol=1e-12)
            for label in [ax.xaxis.label, ax.yaxis.label]:
                expected_size = plots.style.AXIS_SIZE
                check(label.get_fontsize() == expected_size, "Axis font size changed")
                check(label.get_fontfamily() == ["DejaVu Sans"], "Verification text font changed")
            for label in ax.get_legend().get_texts():
                check(label.get_fontsize() >= plots.style.LEGEND_SIZE, "Legend is below 6 pt")
            for spine in ax.spines.values():
                check(spine.get_linewidth() == plots.style.FRAME_WIDTH, "Frame width changed")
            check(ax.get_legend().get_title().get_text() == expected["legend_title"],
                  f"Fig {number}: legend title")
            if number == "2":
                # Check the scientific meaning separately from the original
                # numerical fixture, whose supplied legend was reversed.
                means = [line for line in ax.lines if not line.get_label().startswith("_")]
                check([line.get_label() for line in means] == [
                    r"$\Delta \tilde{S}_{\mathrm{iso}}$", r"$\Delta \tilde{S}$"],
                    "Fig 2C: isotropic/measured surface legend mapping")
                for line, quantity in zip(means, ["volume_scaled_area", "surface_area"]):
                    values = [row[4] for row in rows if row[2] == quantity]
                    np.testing.assert_allclose(line.get_ydata(), values,
                                               rtol=1e-10, atol=1e-12)
                check(" ".join(ax.get_ylabel().split()) == "Normalized surface expansion",
                      "Fig 2C: shared y-axis title")
            check(len(ax.lines) == len(expected["lines"]), f"Fig {number}: line count")
            line_colors = ([color for color, cohort in zip(tol_colors, ["control", "BB94-1uM", "BB94-5uM"])
                            for _ in range(len(plots.COHORTS[cohort]) + 1)]
                           if number == "4" else [line["color"] for line in expected["lines"]])
            for index, (line, wanted) in enumerate(zip(ax.lines, expected["lines"])):
                for values, key in [(line.get_xdata(), "x"), (line.get_ydata(), "y"),
                                    (mpl.colors.to_rgba(line.get_color()), "color")]:
                    target = mpl.colors.to_rgba(line_colors[index]) if key == "color" else wanted[key]
                    np.testing.assert_allclose(values, np.array(target, dtype=float),
                                               rtol=1e-10, atol=1e-12, equal_nan=True,
                                               err_msg=f"Fig {number}, line {index}: {key}")
                is_mean = wanted["linewidth"] == 1
                width = plots.style.MEAN_WIDTH if is_mean else 0.5 if number == "4" else plots.style.DATA_WIDTH
                check(line.get_linewidth() == width, "Line width changed")
                alpha = 0.3 if number == "4" and not is_mean else wanted["alpha"]
                check(line.get_alpha() == alpha, "Line opacity changed")
                line_count += 1
            check(len(ax.collections) == len(expected["bands"]), f"Fig {number}: band count")
            for index, (band, wanted) in enumerate(zip(ax.collections, expected["bands"])):
                check(len(band.get_paths()) == len(wanted["paths"]), "Band path count changed")
                for path, vertices in zip(band.get_paths(), wanted["paths"]):
                    np.testing.assert_allclose(path.vertices, vertices, rtol=1e-10, atol=1e-12)
                facecolors = [mpl.colors.to_rgba(tol_colors[index], .2)] if number == "4" else wanted["facecolors"]
                np.testing.assert_allclose(band.get_facecolors(), facecolors,
                                           rtol=1e-10, atol=1e-12)
                band_count += 1
            stem = args.output_dir / (plots.STEMS[number] + ("-tol-bright" if number == "4" else ""))
            for extension in ["pdf", "svg", "png"]:
                path = stem.with_suffix("." + extension)
                check(path.is_file() and path.stat().st_size > 0, f"Missing export: {path}")
            csv_path = stem.with_name(stem.name + "-summary.csv")
            with csv_path.open(newline="") as stream:
                observed_rows = list(csv.reader(stream))[1:]
            check(len(observed_rows) == len(rows), f"Summary row count changed: {csv_path}")
            for observed, wanted in zip(observed_rows, rows):
                check(observed[:3] == wanted[:3], "Summary identifiers changed")
                np.testing.assert_allclose(np.array(observed[3:], dtype=float), wanted[3:],
                                           rtol=1e-10, atol=1e-12)
                summary_count += 1
            plt.close(fig)
    variant = args.output_dir / (plots.STEMS["4"] + "-tol-bright")
    svg = variant.with_suffix(".svg").read_text()
    svg = re.sub(r"#[0-9a-fA-F]{3,8}\b", lambda m: mpl.colors.to_hex(m.group()).upper(), svg)
    positions = [svg.index(color) for color in tol_colors]
    check(positions == sorted(positions), "Fig 4B: exported Tol color order")
    recolored_curves = len(reference["4"]["lines"])
    recolored_bands = len(reference["4"]["bands"])
    report = {
        "status": "passed", "source_tables_verified": len(sources),
        "source_bytes": sum(s["bytes"] for s in sources),
        "individual_and_mean_curves_verified": line_count, "sd_bands_verified": band_count,
        "summary_rows_verified": summary_count, "plots_verified": len(reference),
        "recolored_curves_verified": recolored_curves,
        "recolored_sd_bands_verified": recolored_bands,
        "fig4b_tol_mapping": {"BB94_um": [0, 1, 5], "matched_z": [.001, .005, .01], "colors": tol_colors},
        "relative_tolerance": 1e-10, "absolute_tolerance": 1e-12,
        "presentation_ylim_overrides": plots.PRESENTATION_YLIMS,
        "fig4b_curve_style": {"individual_line_width_points": 0.5, "individual_opacity": 0.3,
                              "mean_line_width_points": plots.style.MEAN_WIDTH,
                              "mean_opacity": 1, "sd_band_opacity": 0.2},
        "fig2c_legend_mapping": {
            "green": "Delta S-tilde_iso: 100 * ((V/V0)**(2/3) - 1)",
            "orange": "Delta S-tilde: 100 * (A/A0 - 1)",
            "shared_y_axis_label": "Normalized surface expansion",
        },
        "numpy": np.__version__, "matplotlib": mpl.__version__,
        "note": "Every numerical curve and SD polygon agrees with Justin's original scripts. Export existence and CSV values are checked. Vector glyph outlines and exported dimensions are audited separately in the root validation records.",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Verified {len(sources)} mesh tables, {line_count} curves and {band_count} SD bands.")


if __name__ == "__main__":
    main()
