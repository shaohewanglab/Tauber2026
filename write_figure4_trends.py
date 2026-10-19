"""Collect panel-wise Spearman tests; correct the three planned trend tests."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-models", action="store_true")
    args = parser.parse_args()
    reports = json.loads((ROOT/"BB-94-titration/results/combined/Fig4FG-BB94-titration-trends.json").read_text())
    if not args.skip_models:
        reports["Fig4E"] = json.loads((ROOT/"model-morphometry/results/Fig4E-model-bud-count-trend.json").read_text())
    panels = ["Fig4E", "Fig4F", "Fig4G"]
    corrected = multipletests([reports[p]["p_raw"] if p in reports else 1 for p in panels], method="holm")[1]
    for panel, adjusted in zip(panels, corrected):
        if panel in reports:
            reports[panel].update(p_adjusted=float(adjusted), correction="Holm across the three planned Figure 4E/F/G trend tests",
                                  trend_family_size=3)
    report = {"panels": {p: reports[p] for p in panels if p in reports},
              "status": "complete" if len(reports) == 3 else "model trend not regenerated",
              "family": panels, "family_size": 3, "multiplicity": "Holm",
              "separate_from_pairwise_families": True,
              "annotation": "Report trends in the caption; retain the separately specified pairwise comparisons on the plots"}
    out = ROOT/"validation/Fig4EFG-spearman-trends.json"
    out.write_text(json.dumps(report, indent=2)+"\n")
    text = ["Figure 4E/F/G: Spearman permutation trends", "",
            "Two-sided absolute-rho permutation tests; one independent simulation (E) or gland (F/G) per observation.",
            "Holm correction across the three planned trend tests, separately from all pairwise comparison families.", ""]
    for panel in panels:
        if panel in reports:
            r = reports[panel]
            text.append(f"{panel}: n={r['n']}; rho={r['rho']:.6f}; raw P={r['p_raw']:.9g}; Holm P={r['p_adjusted']:.9g}; "
                        f"{r['permutations']} permutations, {r['extreme_permutations']} equally/more extreme; seed={r['seed']}.")
    text.extend(["", "These trend tests assess monotonic association across all thicknesses/doses; they do not identify which individual groups differ from baseline."])
    out.with_suffix(".txt").write_text("\n".join(text)+"\n")
    print("\n".join(text))


if __name__ == "__main__":
    main()
