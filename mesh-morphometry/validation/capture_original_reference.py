#!/usr/bin/env python3
"""Capture numerical curves directly from Justin's supplied plotting scripts.

Supply --scripts PlotScripts_Justin and --paper-data PaperData-v2. Only the
input root and noninteractive output handling are replaced; calculations and
cohort selection run unchanged. No consolidated plotting code is imported.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from unittest.mock import patch

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tauber-mesh-matplotlib"))
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

SCRIPTS = {
    "2": "create_areavolumeratio_plot_experiment.py",
    "4": "create_asphericity_plot_experiment.py",
    "5": "create_asphericity_plot_washout.py",
}


def json_safe(value):
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def capture(scripts, paper_data):
    reference = {}
    sys.path.insert(0, str(scripts))
    for number, name in SCRIPTS.items():
        path = scripts / name
        source, replacements = re.subn(
            r'^rootFolder\s*=\s*"[^"\n]+"',
            lambda _: "rootFolder = " + repr(str(paper_data)),
            path.read_text(), flags=re.MULTILINE,
        )
        if replacements != 1:
            raise ValueError(f"Expected exactly one input root in {path}")
        with mpl.rc_context(), patch.object(plt, "savefig"), patch.object(plt, "show"):
            namespace = {"__name__": "__main__", "__file__": str(path)}
            exec(compile(source, str(path), "exec"), namespace)
            fig, ax = namespace["fig"], namespace["ax"]
            reference[number] = {
                "original_script": "PlotScripts_Justin/" + name,
                "figsize": fig.get_size_inches().tolist(),
                "position": list(ax.get_position().bounds),
                "xlim": list(ax.get_xlim()), "ylim": list(ax.get_ylim()),
                "xlabel": ax.get_xlabel(), "ylabel": ax.get_ylabel(),
                "legend_title": ax.get_legend().get_title().get_text(),
                "lines": [{
                    "x": line.get_xdata().tolist(), "y": line.get_ydata().tolist(),
                    "color": list(mpl.colors.to_rgba(line.get_color())),
                    "alpha": line.get_alpha(), "linewidth": line.get_linewidth(),
                    "label": line.get_label(),
                } for line in ax.lines],
                "bands": [{
                    "paths": [path.vertices.tolist() for path in band.get_paths()],
                    "facecolors": band.get_facecolors().tolist(),
                } for band in ax.collections],
            }
            plt.close(fig)
    return json_safe(reference)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scripts", type=Path, required=True)
    parser.add_argument("--paper-data", type=Path, required=True)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("expected-curves.json"))
    args = parser.parse_args()
    reference = capture(args.scripts.resolve(), args.paper_data.resolve())
    args.output.write_text(json.dumps(reference, separators=(",", ":"), allow_nan=False) + "\n")
    print(f"Captured {sum(len(x['lines']) for x in reference.values())} curves "
          f"and {sum(len(x['bands']) for x in reference.values())} SD bands.")


if __name__ == "__main__":
    main()
