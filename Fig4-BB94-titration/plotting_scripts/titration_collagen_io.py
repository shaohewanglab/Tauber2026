"""Validated line-scan input and gland-level endpoints for collagen IV titration.

Paths in the JSON configuration are relative to that configuration file. Each
CSV describes one dose/cohort and identifies glands with ``file_name`` and
individual scans with ``id``. No scan is treated as an independent replicate.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _nonempty_strings(series: pd.Series, label: str) -> pd.Series:
    cleaned = series.astype("string").str.strip()
    if cleaned.isna().any() or cleaned.eq("").any():
        raise ValueError(f"Missing or empty {label}.")
    return cleaned.astype(str)


def load_collagen_data(config_path):
    """Return a gland endpoint DataFrame and a JSON-serializable audit dict.

    All scans contribute equally within each gland. Each gland then contributes
    equally to the control normalization denominator and subsequent statistics.
    The returned ``value`` is mean gland AUC / mean control gland AUC.
    """
    config_path = Path(config_path).expanduser().resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    if config.get("schema_version") != 1:
        raise ValueError("Collagen configuration must have schema_version 1.")
    inputs = config.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError("Collagen configuration requires a nonempty inputs list.")
    columns = config.get("columns", {})
    gland_column = columns.get("gland", "file_name")
    scan_column = columns.get("scan", "id")
    auc_column = columns.get("auc", "AUC")
    channel_column = columns.get("channel", "channel")
    expected_channel = config.get("expected_channel", 2)
    minimum_scans = config.get("minimum_scans_per_gland", 1)
    if isinstance(minimum_scans, bool) or not isinstance(minimum_scans, int) or minimum_scans < 1:
        raise ValueError("minimum_scans_per_gland must be a positive integer.")

    endpoints = []
    input_audit = []
    control_scan_values = []
    seen_paths = set()
    seen_dose_cohorts = set()
    for entry in inputs:
        if not isinstance(entry, dict):
            raise ValueError("Each collagen input must be an object.")
        source = entry.get("path")
        if not isinstance(source, str) or not source.strip():
            raise ValueError("Every collagen input requires a path.")
        path = (config_path.parent / source).expanduser().resolve()
        if path in seen_paths:
            raise ValueError(f"The same input CSV is listed twice: {path}")
        seen_paths.add(path)
        try:
            dose = float(entry["dose_um"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"Input {path.name} needs a numeric dose_um.") from error
        if not math.isfinite(dose) or dose < 0:
            raise ValueError(f"Invalid dose_um for {path.name}: {dose}")
        cohort = str(entry.get("cohort", "")).strip()
        if not cohort:
            raise ValueError(f"Input {path.name} requires a nonempty cohort.")
        group = str(entry.get("group", "Vehicle" if dose == 0 else f"{dose:g} uM BB-94")).strip()
        if (cohort, dose) in seen_dose_cohorts:
            raise ValueError(f"Duplicate cohort/dose entry ({cohort}, {dose:g}); combine scans in one CSV to prevent duplicate glands.")
        seen_dose_cohorts.add((cohort, dose))
        frame = pd.read_csv(path, dtype={gland_column: "string", scan_column: "string"})
        required = {gland_column, scan_column, auc_column, channel_column}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
        if frame.empty:
            raise ValueError(f"No line scans in {path.name}.")
        frame[gland_column] = _nonempty_strings(frame[gland_column], f"gland IDs in {path.name}")
        frame[scan_column] = _nonempty_strings(frame[scan_column], f"scan IDs in {path.name}")
        if frame.duplicated([gland_column, scan_column]).any():
            raise ValueError(f"Duplicate gland/scan IDs in {path.name}.")
        auc = pd.to_numeric(frame[auc_column], errors="coerce")
        if not np.isfinite(auc.to_numpy(dtype=float)).all() or auc.lt(0).any():
            raise ValueError(f"AUC must be finite and nonnegative in {path.name}.")
        frame[auc_column] = auc.astype(float)
        channels = pd.to_numeric(frame[channel_column], errors="coerce")
        if not channels.eq(expected_channel).all():
            raise ValueError(f"All scans must be channel {expected_channel} in {path.name}.")
        gland_means = frame.groupby(gland_column, sort=True)[auc_column].agg(["mean", "size"])
        if gland_means["size"].lt(minimum_scans).any():
            failures = gland_means.index[gland_means["size"].lt(minimum_scans)].tolist()
            raise ValueError(f"{path.name}: glands {failures} have fewer than {minimum_scans} scans.")
        expected_glands = entry.get("expected_glands")
        if expected_glands is not None:
            if isinstance(expected_glands, bool) or not isinstance(expected_glands, int) or expected_glands < 1:
                raise ValueError(f"expected_glands must be a positive integer for {path.name}.")
            if len(gland_means) != expected_glands:
                raise ValueError(f"{path.name}: expected {expected_glands} glands, found {len(gland_means)}.")
        for source_gland, row in gland_means.iterrows():
            gland_id = f"{cohort}/dose-{dose:g}/gland-{source_gland}"
            endpoints.append({
                "dose_um": dose, "cohort": cohort, "gland_id": gland_id,
                "source_gland_id": str(source_gland), "group": group,
                "raw_auc": float(row["mean"]), "n_scans": int(row["size"]),
                "source_file": str(path),
            })
        if dose == 0:
            control_scan_values.extend(frame[auc_column].tolist())
        input_audit.append({
            "path": str(path), "sha256": _sha256(path), "cohort": cohort,
            "dose_um": dose, "group": group, "n_scans": len(frame),
            "n_glands": len(gland_means),
            "scans_per_gland": {str(key): int(value) for key, value in gland_means["size"].items()},
        })

    result = pd.DataFrame(endpoints).sort_values(["dose_um", "cohort", "source_gland_id"]).reset_index(drop=True)
    if result["gland_id"].duplicated().any():
        raise ValueError("Duplicate gland identities across collagen inputs.")
    controls = result.loc[result["dose_um"].eq(0), "raw_auc"]
    if controls.empty:
        raise ValueError("A vehicle-control input with dose_um=0 is required.")
    denominator = float(controls.mean())
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("The mean control gland AUC must be positive and finite.")
    result["value"] = result["raw_auc"] / denominator
    weighted_denominator = float(np.mean(control_scan_values))
    audit = {
        "analysis": "collagen_iv_intensity", "statistical_unit": "gland",
        "endpoint": "mean raw-image collagen IV line-scan AUC per gland, normalized to mean vehicle gland AUC",
        "config_path": str(config_path), "config_sha256": _sha256(config_path),
        "inputs": input_audit, "n_glands": len(result),
        "n_scans": int(result["n_scans"].sum()),
        "minimum_scans_per_gland": int(result["n_scans"].min()),
        "maximum_scans_per_gland": int(result["n_scans"].max()),
        "excluded_glands": [], "excluded_scans": [],
        "normalization": {
            "method": "Equal-weight mean of vehicle-control gland means",
            "control_glands": len(controls), "control_gland_mean_auc": denominator,
            "notebook_pooled_control_scan_mean_auc": weighted_denominator,
            "new_value_divided_by_notebook_value": weighted_denominator / denominator,
            "note": "The reference notebook and earlier gland-level correction used the pooled control-scan mean. Equal gland weighting makes the mean of plotted control glands exactly 1. This common scale factor does not change test statistics or p-values.",
        },
    }
    return result, audit
