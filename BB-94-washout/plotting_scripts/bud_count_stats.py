"""Validated gland-level interval preparation, extracted from the washout workflow.

Only matching and ratio calculation are retained. Inferential comparisons for
this package are Mann-Whitney U with Holm correction, in the plotting module.
"""
import json
import math
import numpy as np
import pandas as pd


def _require_columns(frame, columns, label):
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError(f"{label} is empty.")


def _number_label(value):
    return f"{float(value):g}"


def _positive_numeric(series, label, integer=False):
    values = pd.to_numeric(series, errors="coerce")
    if not np.isfinite(values.to_numpy(dtype=float)).all() or (values <= 0).any():
        raise ValueError(f"{label} must contain finite, positive values; no pseudocount is added.")
    if integer and not np.equal(values, np.floor(values)).all():
        raise ValueError(f"{label} must contain integer bud counts.")
    return values


def prepare_intervals(decoded, times=(0, 22, 46), group_order=("ctrl", "BB94")):
    """Match each gland across three times and return its two consecutive ratios.

    Input columns: experiment, group, gland, time_h, bud_number. An optional
    subject_id is preserved after checking that it identifies exactly one gland.
    All glands must have exactly one observation at every requested time. Dates,
    row order, and equal gland numbers in different treatments never create pairs.
    """
    required = ["experiment", "group", "gland", "time_h", "bud_number"]
    _require_columns(decoded, required, "Decoded counts")
    if len(times) != 3:
        raise ValueError("Exactly three observation times are required.")
    try:
        times = tuple(float(t) for t in times)
    except (TypeError, ValueError) as exc:
        raise ValueError("Observation times must be numeric.") from exc
    if not all(math.isfinite(t) for t in times) or not times[0] < times[1] < times[2]:
        raise ValueError("Observation times must be finite and strictly increasing.")
    group_order = tuple(str(g) for g in group_order)
    if len(group_order) != 2 or len(set(group_order)) != 2:
        raise ValueError("Exactly two distinct treatment groups are required.")

    frame = decoded.copy()
    identity = ["experiment", "group", "gland"]
    for column in identity:
        if frame[column].isna().any():
            raise ValueError(f"{column} contains a missing identifier.")
        frame[column] = frame[column].astype(str).str.strip()
        if frame[column].eq("").any():
            raise ValueError(f"{column} contains an empty identifier.")
    if set(frame["group"]) != set(group_order):
        raise ValueError(f"Expected exactly the treatment groups {group_order}.")
    frame["time_h"] = pd.to_numeric(frame["time_h"], errors="coerce")
    if not np.isfinite(frame["time_h"].to_numpy(dtype=float)).all():
        raise ValueError("time_h must contain finite numeric values.")
    if set(frame["time_h"]) != set(times):
        raise ValueError(f"Expected exactly the observation times {times}.")
    frame["bud_number"] = _positive_numeric(frame["bud_number"], "bud_number", integer=True)
    if frame.duplicated(identity + ["time_h"]).any():
        raise ValueError("Duplicate count for the same experiment, treatment, gland, and time.")

    if "subject_id" in frame:
        if frame["subject_id"].isna().any() or frame["subject_id"].astype(str).str.strip().eq("").any():
            raise ValueError("subject_id contains a missing or empty identifier.")
        frame["subject_id"] = frame["subject_id"].astype(str)
        mappings = frame[identity + ["subject_id"]].drop_duplicates()
        if mappings.duplicated(identity).any() or mappings.duplicated("subject_id").any():
            raise ValueError("subject_id must map one-to-one to experiment, treatment, and gland.")
    else:
        frame["subject_id"] = [
            json.dumps(list(key), ensure_ascii=False, separators=(",", ":"))
            for key in frame[identity].itertuples(index=False, name=None)
        ]

    records = []
    for key, observations in frame.groupby(identity, sort=True, observed=True):
        if set(observations["time_h"]) != set(times):
            missing = sorted(set(times) - set(observations["time_h"]))
            raise ValueError(f"Gland {key} is missing observation time(s) {missing}; incomplete glands are not dropped.")
        counts = observations.set_index("time_h")["bud_number"]
        for start, end in zip(times[:-1], times[1:]):
            bud_start, bud_end = int(counts.loc[start]), int(counts.loc[end])
            ratio = bud_end / bud_start
            records.append({
                "experiment": key[0], "group": key[1], "gland": key[2],
                "subject_id": observations["subject_id"].iloc[0],
                "start_h": start, "end_h": end,
                "interval": f"{_number_label(end)}/{_number_label(start)}",
                "bud_start": bud_start, "bud_end": bud_end,
                "ratio": ratio, "log2_ratio": float(np.log2(ratio)),
            })
    result = pd.DataFrame.from_records(records)
    result["_group_order"] = result["group"].map(dict(zip(group_order, range(2))))
    result = result.sort_values(["start_h", "_group_order", "experiment", "gland"]).drop(columns="_group_order").reset_index(drop=True)
    result.attrs["group_order"] = group_order
    result.attrs["times"] = times
    return result
