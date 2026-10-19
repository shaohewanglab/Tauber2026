#!/usr/bin/env python3
"""Reproduce all manuscript experimental/model plots from bundled inputs."""
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

    for results in [ROOT/"BB-94-titration/results", ROOT/"BB-94-washout/results",
                    ROOT/"model-morphometry/results"]:
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
    parser.add_argument("--mesh-text-rendering", choices=["latex", "mathtext"], default="latex",
                        help="XeLaTeX retains Justin's Arev math with Arial text; mathtext needs no TeX.")
    parser.add_argument("--skip-models", action="store_true",
                        help="Rebuild experimental panels only; model plots require XeLaTeX and Arial.")
    args = parser.parse_args()
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir())/"tauber-modeling-matplotlib"))
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    fig4 = ROOT/"BB-94-titration/plotting_scripts"
    fig5 = ROOT/"BB-94-washout/plotting_scripts"
    commands = [
        [fig4/"plot_bud_count_titration.py", "--output-dir", ROOT/"BB-94-titration/results/bud-count", "--font", args.font, "--no-plot", "--overwrite"],
        [fig4/"plot_collagen_iv_titration.py", "--output-dir", ROOT/"BB-94-titration/results/collagen-IV", "--font", args.font, "--no-plot", "--overwrite"],
        [ROOT/"reproduce_combined_titration.py", "--font", args.font],
        [fig5/"prepare_washout_data.py"],
        [fig5/"plot_combined_log2_ratios_uniform.py", "--font", args.font],
        [fig5/"plot_combined_log2_ratios_symbols_no_legend.py", "--font", args.font],
        [ROOT/"mesh-morphometry/reproduce_mesh_plots.py",
         "--text-rendering", args.mesh_text_rendering, "--font", args.font],
    ]
    if not args.skip_models:
        commands.append([ROOT/"model-morphometry/reproduce_model_plots.py"])
    for command in commands:
        print("Running " + Path(command[0]).name, flush=True)
        subprocess.run([sys.executable, *map(str, command)], cwd=ROOT, env=env, check=True)
        # Normalize each producer's paths before later plotters hash its tables.
        portable_provenance()
    trends = [sys.executable, str(ROOT/"write_figure4_trends.py")]
    if args.skip_models:
        trends.append("--skip-models")
    subprocess.run(trends, cwd=ROOT, env=env, check=True)
    verify = [sys.executable, str(ROOT/"verify_reproduction.py")]
    if args.skip_models:
        verify.append("--skip-models")
    subprocess.run(verify, cwd=ROOT, env=env, check=True)
    print("Reproduction complete. See the results folders for PDF, SVG, PNG and statistical reports.")


if __name__ == "__main__":
    main()
