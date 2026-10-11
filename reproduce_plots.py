#!/usr/bin/env python3
"""Reproduce the current Fig. 4 and Fig. 5 plots from the included source tables."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def portable_provenance():
    """Store paths within this package relatively in generated audit files."""
    prefix = str(ROOT) + os.sep

    def clean(value):
        if isinstance(value, str):
            return value.replace(prefix, "")
        if isinstance(value, list):
            return [clean(item) for item in value]
        if isinstance(value, dict):
            return {key: clean(item) for key, item in value.items()}
        return value

    for results in [ROOT/"Fig4-BB94-titration/results", ROOT/"Fig5-washout/results"]:
        for path in results.rglob("*"):
            if path.suffix == ".json":
                record = json.loads(path.read_text())
                path.write_text(json.dumps(clean(record), indent=2, allow_nan=False)+"\n")
            elif path.suffix in {".csv", ".txt", ".html"}:
                original = path.read_text()
                changed = original.replace(prefix, "")
                if changed != original:
                    path.write_text(changed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--font", default="Arial", help="Use Arial for the manuscript appearance; e.g. DejaVu Sans for systems without Arial.")
    args = parser.parse_args()
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir())/"tauber-modeling-matplotlib"))
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    fig4 = ROOT/"Fig4-BB94-titration/plotting_scripts"
    fig5 = ROOT/"Fig5-washout/plotting_scripts"
    commands = [
        [fig4/"plot_bud_count_titration.py", "--output-dir", ROOT/"Fig4-BB94-titration/results/bud-count", "--font", args.font, "--overwrite"],
        [fig4/"plot_collagen_iv_titration.py", "--output-dir", ROOT/"Fig4-BB94-titration/results/collagen-IV", "--font", args.font, "--overwrite"],
        [fig5/"prepare_washout_data.py"],
        [fig5/"plot_combined_log2_ratios_uniform.py", "--font", args.font],
        [fig5/"plot_combined_log2_ratios_symbols_no_legend.py", "--font", args.font],
    ]
    for command in commands:
        print("Running " + Path(command[0]).name, flush=True)
        subprocess.run([sys.executable, *map(str, command)], cwd=ROOT, env=env, check=True)
    portable_provenance()
    subprocess.run([sys.executable, str(ROOT/"verify_reproduction.py")], cwd=ROOT, env=env, check=True)
    print("Reproduction complete. See the results folders for PDF, SVG, PNG and statistical reports.")


if __name__ == "__main__":
    main()
