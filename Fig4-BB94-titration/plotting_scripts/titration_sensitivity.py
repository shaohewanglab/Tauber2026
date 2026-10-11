"""Cohort-adjusted sensitivity analysis for independent gland-level responses.

This supplements the pooled all-pairs analysis; its results do not determine the
asterisks in the primary plot. The outcome is one log2 count ratio per gland.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from patsy import build_design_matrices
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests


CONTRAST_COLUMNS = [
    "dose_um", "control_dose_um", "mean_difference", "standard_error",
    "t_statistic", "df", "p_unadjusted", "p_adjusted", "ci95_low",
    "ci95_high", "reject",
]


def bud_cohort_sensitivity(data: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Fit dose + cohort OLS and compare every dose with the vehicle control.

    HC3 covariance and residual-df t inference are explicit. Holm correction
    covers the dose-versus-control contrasts only. The confidence intervals are
    unadjusted 95% intervals, not simultaneous intervals. Design-vector
    contrasts avoid depending on coefficient names or categorical reference
    levels. Each gland must contribute exactly one outcome.
    """
    required = {"dose_um", "value", "cohort", "gland_id"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing sensitivity-analysis columns: {sorted(missing)}")
    frame = data.loc[:, sorted(required)].copy()
    if frame.empty or frame.isna().any().any():
        raise ValueError("Sensitivity analysis requires complete, nonempty gland data.")
    if frame["gland_id"].duplicated().any():
        raise ValueError("Sensitivity analysis requires one response per gland_id.")
    for name in ["dose_um", "value"]:
        frame[name] = pd.to_numeric(frame[name], errors="raise")
        if not np.isfinite(frame[name]).all():
            raise ValueError(f"{name} must contain finite numbers.")
    frame["cohort"] = frame["cohort"].astype(str)
    doses = sorted(frame["dose_um"].unique().tolist())
    cohorts = sorted(frame["cohort"].unique().tolist())
    metadata = {
        "available": False,
        "method": "OLS: per-gland log2 count ratio ~ categorical dose + categorical cohort",
        "covariance": "HC3",
        "use_t": True,
        "n_glands": len(frame),
        "n_cohorts": len(cohorts),
        "cohort_levels": cohorts,
        "control_dose_um": 0.0,
        "correction": "Holm, across dose-versus-control contrasts",
        "confidence_intervals": "Unadjusted 95% t intervals for dose minus control",
        "purpose": "Sensitivity analysis only; does not determine primary plot asterisks",
    }
    if 0.0 not in doses or len(doses) < 2:
        raise ValueError("Sensitivity analysis requires vehicle (0 µM) and a treatment.")
    if len(cohorts) < 2:
        metadata["reason"] = "At least two cohorts are needed for cohort adjustment."
        return pd.DataFrame(columns=CONTRAST_COLUMNS), metadata

    model = ols("value ~ C(dose_um) + C(cohort)", data=frame).fit(
        cov_type="HC3", use_t=True
    )
    design = model.model.exog
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("Dose and cohort are confounded; cohort-adjusted contrasts are not identifiable.")
    if model.df_resid <= 0:
        raise ValueError("Insufficient residual degrees of freedom for cohort adjustment.")
    if not np.isfinite(model.cov_params()).all().all():
        raise ValueError("HC3 covariance is not finite; check replication within dose/cohort groups.")

    def average_design(dose: float) -> np.ndarray:
        prediction_frame = pd.DataFrame({"dose_um": dose, "cohort": cohorts})
        matrix = build_design_matrices(
            [model.model.data.design_info], prediction_frame
        )[0]
        return np.asarray(matrix).mean(axis=0)

    control = average_design(0.0)
    rows = []
    for dose in doses:
        if dose == 0.0:
            continue
        contrast = model.t_test(average_design(dose) - control, use_t=True)
        interval = np.asarray(contrast.conf_int(alpha=0.05)).reshape(-1)
        scalar = lambda value: float(np.asarray(value).reshape(-1)[0])
        rows.append({
            "dose_um": float(dose),
            "control_dose_um": 0.0,
            "mean_difference": scalar(contrast.effect),
            "standard_error": scalar(contrast.sd),
            "t_statistic": scalar(contrast.tvalue),
            "df": float(model.df_resid),
            "p_unadjusted": scalar(contrast.pvalue),
            "ci95_low": float(interval[0]),
            "ci95_high": float(interval[1]),
        })
    result = pd.DataFrame(rows)
    reject, adjusted, _, _ = multipletests(result["p_unadjusted"], method="holm")
    result["p_adjusted"] = adjusted
    result["reject"] = reject
    metadata.update({
        "available": True,
        "n_comparisons": len(result),
        "residual_df": float(model.df_resid),
        "adjustment": "Common additive cohort effect; no dose-by-cohort interaction fitted",
    })
    return result.loc[:, CONTRAST_COLUMNS], metadata
