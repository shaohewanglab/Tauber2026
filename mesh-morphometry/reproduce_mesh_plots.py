#!/usr/bin/env python3
"""Reproduce Justin's three cleaned-mesh time-course plots from 28 small tables.

Adapted from Justin's create_areavolumeratio_plot_experiment.py,
create_asphericity_plot_experiment.py, create_asphericity_plot_washout.py,
and the PNAS settings in PlotLibrary.py. Numerical operations, cohort order,
are retained; Figure 4 uses the selected Tol Bright colors, the reversed
Figure 2C legend labels are corrected, and presentation uses the shared style.
"""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tauber-mesh-matplotlib"))

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_style as style

ROOT = Path(__file__).resolve().parent
COHORTS = {
    "control": ["241104-SMG2", "241104-SMG4", "241104-SMG5", "241105-SMG1",
                "241105-SMG2", "241105-SMG3", "241105-SMG4_SG1"],
    "BB94-1uM": [f"241114-SMG{i}-BB94-1uM" for i in range(1, 6)],
    "BB94-5uM": [f"250401-SMG{i}-BB94-5uM" for i in range(1, 9)],
    "washout": [f"250916-SMG{i}-BB94-washout" for i in range(1, 9)],
}
STEMS = {
    "2": "Fig2C-area_volume_evolution_control_cleanmesh_plus",
    "4": "Fig4B-asphericitychange_evolution_BB94_cleanmesh_plus",
    "5": "Fig5F-asphericitychange_evolution_washout_cleanmesh_plus",
}
# Display limits are independent of the frozen original numerical curves.
PRESENTATION_YLIMS = {"4": (-10, 110)}
def configure_style(text_rendering="latex", font="Arial"):
    """Keep Justin's palette and math font, with common Arial typography."""
    plt.style.use("seaborn-v0_8-colorblind")
    mpl.rcParams.update(style.rc(font="Arial" if text_rendering == "latex" else font,
                                arev_math=text_rendering == "latex"))
    mpl.rcParams.update({
        "xtick.major.pad": 1, "xtick.minor.pad": 1,
        "ytick.major.pad": 1, "ytick.minor.pad": 1,
        "axes.labelpad": style.MESH_AXIS_LABEL_PAD,
        "legend.labelspacing": 0.3, "axes.axisbelow": False,
    })


def load_cohort(name):
    result = []
    for sample in COHORTS[name]:
        path = ROOT / "source_data" / f"volume_area_{sample}_cp_model0907_v2.txt"
        data = np.genfromtxt(path, skip_header=1)
        if data.ndim != 2 or data.shape[1] != 4 or not np.isfinite(data[:, :2]).all():
            raise ValueError(f"Invalid frame/time/volume/area table: {path}")
        if not np.all(np.diff(data[:, 1]) > 0) or data[0, 1] != 0:
            raise ValueError(f"Invalid time axis: {path}")
        measures = data[:, 2:]
        if not np.isfinite(measures[0]).all() or np.any(measures[np.isfinite(measures)] <= 0):
            raise ValueError(f"Nonpositive or missing baseline measurement: {path}")
        # The original NaN volume is retained and excluded only by nanmean/nanstd.
        result.append((sample, data[:, 1] / 60, data[:, 2], data[:, 3]))
    return result


def asphericity(area, volume):
    sphericity = (np.pi ** (1 / 3) * (6 * volume) ** (2 / 3)) / area
    return 1 / sphericity - 1


def summarize(times, curves):
    # Match Justin: pad short series with NaNs and average by row, without
    # interpolation. A few timestamps differ by <=5 min. Use the first longest
    # time axis, preserving Python's max(..., key=len) tie behavior.
    padded = np.full((len(curves), max(map(len, curves))), np.nan)
    for index, curve in enumerate(curves):
        padded[index, :len(curve)] = curve
    return (max(times, key=len), np.nanmean(padded, axis=0),
            np.nanstd(padded, axis=0), np.isfinite(padded).sum(axis=0))


def add_mean(ax, figure, cohort, quantity, times, curves, color, label, rows,
             scale_before_mean=True):
    # Fig. 2 multiplies by 100 after padding; asphericity multiplies each curve
    # first. Keep that operation order for agreement with the original scripts.
    if not scale_before_mean:
        curves = [100 * curve for curve in curves]
    time, mean, sd, n = summarize(times, curves)
    ax.plot(time, mean, color=color, lw=style.MEAN_WIDTH, label=label)
    ax.fill_between(time, mean - sd, mean + sd, color=color, alpha=0.2, lw=0)
    for t, m, s, count in zip(time, mean, sd, n):
        rows.append([figure, cohort, quantity, t, m, s, int(count)])


def format_axes(ax, position, ylim, ymajor, yminor, ylabel, legend_position):
    ax.tick_params("both", left=True, right=False, bottom=True, top=False,
                   direction="out", width=style.FRAME_WIDTH, length=3, which="major")
    ax.tick_params("both", left=True, right=False, bottom=True, top=False,
                   direction="out", width=style.FRAME_WIDTH, length=1.5, which="minor")
    for spine in ax.spines.values():
        spine.set_linewidth(style.FRAME_WIDTH)
    ax.set_xlabel(r"$t$ [h]")
    ax.set_ylabel(ylabel)
    legend = ax.legend(frameon=False, loc=legend_position, handletextpad=1,
                       handlelength=0.6)
    legend.get_frame().set_facecolor("none")
    ax.set_xlim(0, 16)
    ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(ymajor))
    ax.yaxis.set_minor_locator(MultipleLocator(yminor))
    ax.set_position(position)
    return legend


def build_figure(figure, text_rendering="latex", layout=True, font="Arial", palette="tol-bright"):
    """Return the figure and its full-precision mean/SD/sample-size rows."""
    configure_style(text_rendering, font)
    rows = []
    fig, ax = plt.subplots(figsize=style.dimensions(f"fig{figure}-mesh"))
    subscript = r"\text" if text_rendering == "latex" else r"\mathrm"
    if figure == "2":
        times, volume_curves, area_curves = [], [], []
        for _, time, volume, area in load_cohort("control"):
            dv = (volume / volume[0]) ** (2 / 3) - 1
            da = area / area[0] - 1
            times.append(time)
            volume_curves.append(dv)
            area_curves.append(da)
            ax.plot(time, 100 * da, color="C02", lw=style.DATA_WIDTH, alpha=0.5)
            ax.plot(time, 100 * dv, color="C01", lw=style.DATA_WIDTH, alpha=0.5)
        # Isotropic surface expansion follows V**(2/3); measured expansion
        # follows A. The supplied script had these two legend labels reversed.
        add_mean(ax, figure, "control", "volume_scaled_area", times, volume_curves,
                 "C01", rf"$\Delta \tilde{{S}}_{{{subscript}{{iso}}}}$", rows, scale_before_mean=False)
        add_mean(ax, figure, "control", "surface_area", times, area_curves,
                 "C02", r"$\Delta \tilde{S}$", rows, scale_before_mean=False)
        format_axes(ax, [0.15, 0.20, 0.80, 0.70], [-5, 75], 25, 5,
                    "Normalized surface\nexpansion", 0)
    else:
        if figure == "4":
            labels = ["0 µM", "1 µM", "5 µM"]
            colors = style.thickness_colors(palette,
                        style.palettes()["thickness"]["fig4b"]["matched_z"])
            groups = list(zip(["control", "BB94-1uM", "BB94-5uM"],
                              colors, labels))
        else:
            groups = [("washout", style.palettes()["washout"]["BB94"], "BB-94 washout"),
                      ("control", "darkgrey", "Control")]
        for cohort, color, label in groups:
            times, curves = [], []
            for _, time, volume, area in load_cohort(cohort):
                chi = asphericity(area, volume)
                y = 100 * (chi - chi[0])
                times.append(time)
                curves.append(y)
                ax.plot(time, y, color=color,
                        lw=0.5 if figure == "4" else style.DATA_WIDTH,
                        alpha=0.3 if figure == "4" else 0.5)
            add_mean(ax, figure, cohort, "asphericity", times, curves, color, label, rows)
        if figure == "4":
            legend = format_axes(ax, [0.17, 0.27, 0.78, 0.67], PRESENTATION_YLIMS["4"], 50, 10,
                                 rf"$\Delta\chi_{{{subscript}{{3D}}}}$ [%]", 2)
            legend.set_title("[BB-94]:")
            legend.borderaxespad = 0
        else:
            format_axes(ax, [0.10, 0.20, 0.85, 0.70], [-10, 200], 50, 10,
                        rf"$\Delta\chi_{{{subscript}{{3D}}}}$ [%]", 2)
    if layout:
        style.compact_layout(fig, ax, f"fig{figure}-mesh",
                             arev_math=text_rendering == "latex",
                             vector_mathtext=text_rendering == "mathtext")
    return fig, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", choices=["all", *STEMS], default="all")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--palette", choices=["tol-bright"], default="tol-bright",
                        help="Figure 4B uses Paul Tol Bright. Figures 2C and 5F retain their colors.")
    parser.add_argument("--font", default="Arial", help="Text font for the mathtext fallback; XeLaTeX output uses Arial.")
    parser.add_argument("--text-rendering", choices=["latex", "mathtext"], default="latex",
                        help="XeLaTeX combines Arial text with Justin's Arev math; mathtext needs no TeX.")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for number in (STEMS if args.figure == "all" else [args.figure]):
        with mpl.rc_context():
            fig, rows = build_figure(number, args.text_rendering, font=args.font, palette=args.palette)
            suffix = "-tol-bright" if number == "4" else ""
            stem = args.output_dir / (STEMS[number] + suffix)
            if args.text_rendering == "latex":
                style.save_arev_figure(fig, stem)
            else:
                style.save_native_figure(fig, stem, transparent=True)
            with stem.with_name(stem.name + "-summary.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(["figure", "cohort", "quantity", "time_h", "mean_percent",
                                 "sd_percent", "n"])
                writer.writerows(rows)
            plt.close(fig)
            print(f"Reproduced {stem.name}", flush=True)


if __name__ == "__main__":
    main()
