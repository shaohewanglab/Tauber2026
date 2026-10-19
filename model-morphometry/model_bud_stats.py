"""Independent-simulation Dunnett comparisons against the thinnest z."""
from __future__ import annotations

import numpy as np
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dunnett_statistics import dunnett_vs_control

MIN_REPLICATES = 2  # Same replication requirement as experimental 4F/G.


def significance(p):
    if not np.isfinite(p) or not 0 <= p <= 1:
        raise ValueError("Significance requires a finite p-value between 0 and 1.")
    for threshold, label in [(1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (.05, "*")]:
        if p < threshold:
            return label
    return "ns"


def analyze_bud_counts(data, z_order):
    """Compare mean bud counts using the full planned Dunnett family.

    Incomplete families remain n/a until every group has enough independent
    simulations; no missing observations or variances are imputed.
    """
    required = {"z", "replicate", "N_buds"}
    if not required <= set(data):
        raise ValueError(f"Missing columns: {sorted(required-set(data))}")
    if data[list(required)].isna().any().any() or not np.isfinite(data[list(required)].to_numpy(float)).all():
        raise ValueError("Every observed simulation needs finite z, identity and bud count.")
    if data.duplicated(["z", "replicate"]).any():
        raise ValueError("Each independent simulation must contribute once per z.")
    if (data.N_buds < 1).any() or (data.N_buds != np.floor(data.N_buds)).any():
        raise ValueError("Bud counts must be integers at least 1; an unbudded shape counts as 1.")
    order = sorted(z_order)
    if len(order) < 2 or len(set(order)) != len(order) or set(data.z) != set(order):
        raise ValueError("Each planned thickness must have observed simulations, without extra groups.")
    baseline = data.loc[data.z == order[0], "N_buds"].to_numpy(float)
    samples = [data.loc[data.z == z, "N_buds"].to_numpy(float) for z in order[1:]]
    complete = min(map(len, [baseline, *samples])) >= MIN_REPLICATES
    test = dunnett_vs_control(samples, baseline) if complete else None
    rows = []
    for i, (z, sample) in enumerate(zip(order[1:], samples)):
        row = {"baseline_z": order[0], "comparison_z": z,
               "n_baseline": len(baseline), "n_comparison": len(sample),
               "mean_difference_comparison_minus_baseline": float(sample.mean()-baseline.mean()),
               "t_comparison": np.nan, "df": np.nan, "pooled_variance": np.nan,
               "p_raw": np.nan, "p_adjusted": np.nan,
               "p_method": "not_tested", "family_size": len(order)-1,
               "status": "incomplete_family", "annotation": "n/a"}
        if test is not None:
            row.update(t_comparison=float(test["statistic"][i]), df=test["df"],
                       pooled_variance=test["pooled_variance"],
                       p_raw=float(test["p_raw"][i]), p_adjusted=float(test["p_adjusted"][i]),
                       p_method="two-sided single-step Dunnett", status="tested",
                       annotation=significance(test["p_adjusted"][i]))
        rows.append(row)
    return pd.DataFrame(rows)
