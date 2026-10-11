"""Decode blinded titration counts and construct one paired endpoint per gland.

The JSON input manifest supplies dose assignments and acquisition days explicitly.
Rows are joined on randomized filenames; days are paired on cohort/group/gland
identity, never on sort order. This module does not modify any input files.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_table(path: Path, required: tuple[str, ...]) -> pd.DataFrame:
    text = path.read_text(encoding="utf-8-sig")
    if not text.strip():
        raise ValueError(f"Empty input file: {path}")
    try:
        delimiter = csv.Sniffer().sniff(text[:8192], delimiters="\t,").delimiter
    except csv.Error as exc:
        raise ValueError(f"Expected a tab- or comma-delimited table: {path}") from exc
    df = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False)
    df.columns = df.columns.str.strip()
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    df = df.loc[:, list(required)].copy()
    for column in required:
        df[column] = df[column].str.strip()
        if df[column].eq("").any():
            raise ValueError(f"{path}: empty {column}")
    return df


def _decode_pair(record_path: Path, count_path: Path) -> pd.DataFrame:
    keys = _read_table(record_path, ("original_file_name", "converted_file_name"))
    counts = _read_table(count_path, ("file_name", "bud_number"))
    for frame, columns, path in (
        (keys, ("original_file_name", "converted_file_name"), record_path),
        (counts, ("file_name",), count_path),
    ):
        for column in columns:
            if frame[column].duplicated().any():
                raise ValueError(f"{path}: duplicate {column}; pairing must be one-to-one")
    expected = set(keys.converted_file_name)
    observed = set(counts.file_name)
    if expected != observed:
        raise ValueError(
            f"Filename-record/count mismatch for {count_path}: "
            f"missing counts={sorted(expected - observed)}, "
            f"unmatched counts={sorted(observed - expected)}"
        )
    counts["bud_number"] = pd.to_numeric(counts.bud_number, errors="raise")
    values = counts.bud_number.to_numpy(dtype=float)
    if not (np.isfinite(values).all() and (values >= 0).all() and (values == np.floor(values)).all()):
        raise ValueError(f"{count_path}: bud_number must contain finite, nonnegative integers")
    counts["bud_number"] = counts.bud_number.astype(int)
    return keys.merge(counts, left_on="converted_file_name", right_on="file_name", validate="one_to_one")


def _parse_name(filename: str, cohort: dict, source: dict) -> dict:
    kind = cohort["filename_format"]
    name = Path(filename).name
    if name != filename:
        raise ValueError(f"Expected filenames without directory components: {filename}")
    if kind == "treatment_gland":
        match = re.fullmatch(
            r"(?P<group>\d+)-(?:(?:bb(?:-?94)?)-)?"
            r"(?P<treatment>dmso|\d+(?:[p.]\d+)?[unµμ]m)_"
            r"(?P<gland>\d+)(?P<placeholder>-space-holder)?\.[^.]+",
            name, re.IGNORECASE,
        )
        if not match:
            raise ValueError(f"Cannot parse treatment/gland filename: {filename}")
        if "day" not in source:
            raise ValueError("Each treatment_gland source must specify its day")
        day = source["day"]
    elif kind == "date_group_gland":
        match = re.fullmatch(
            r"(?P<date>\d{6})-G(?P<group>\d+)_(?P<gland>\d+)"
            r"(?P<placeholder>-space-holder)?\.[^.]+", name, re.IGNORECASE,
        )
        if not match:
            raise ValueError(f"Cannot parse date/group/gland filename: {filename}")
        date = match.group("date")
        mapping = cohort.get("date_to_day", {})
        if date not in mapping:
            raise ValueError(f"{filename}: date {date} is absent from date_to_day")
        day = mapping[date]
        if "day" in source and source["day"] != day:
            raise ValueError(f"{filename}: source day disagrees with date_to_day")
    else:
        raise ValueError(f"Unknown filename_format {kind!r}")
    if isinstance(day, bool) or not isinstance(day, int):
        raise ValueError(f"{filename}: acquisition day must be an integer")
    group = str(int(match.group("group")))
    dose_map = {str(int(key)): value for key, value in cohort["dose_um_by_group"].items()}
    if group not in dose_map:
        raise ValueError(f"{filename}: group {group} is absent from dose_um_by_group")
    dose = float(dose_map[group])
    if not math.isfinite(dose) or dose < 0:
        raise ValueError(f"Invalid dose for group {group}: {dose}")
    if kind == "treatment_gland":
        treatment = match.group("treatment").lower().replace("p", ".")
        named_dose = 0.0 if treatment == "dmso" else float(treatment[:-2]) / (1000 if treatment.endswith("nm") else 1)
        if not math.isclose(dose, named_dose, abs_tol=1e-12):
            raise ValueError(f"{filename}: filename dose {named_dose} disagrees with manifest dose {dose}")
    gland = str(int(match.group("gland")))
    return {
        "group_id": group,
        "gland_number": gland,
        "gland_id": f"{cohort['cohort']}:G{int(group):02d}:{int(gland):04d}",
        "dose_um": dose,
        "day": day,
        "placeholder": bool(match.group("placeholder")),
    }


def load_bud_data(config_path: str | Path) -> tuple[pd.DataFrame, dict]:
    """Return gland log2(day1/day0) endpoints and JSON-serializable provenance.

    The endpoint days may be changed in the manifest. Filename-marked placeholders
    are excluded only when explicitly enabled for that cohort; all exclusions and
    non-endpoint days are recorded. Counts <=1 are never silently filtered.
    """
    config_path = Path(config_path).expanduser().resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != 1:
        raise ValueError("Bud input manifest requires schema_version 1")
    baseline_day = config.get("baseline_day", 0)
    final_day = config.get("final_day", 1)
    if any(isinstance(day, bool) or not isinstance(day, int) for day in (baseline_day, final_day)) or baseline_day >= final_day:
        raise ValueError("baseline_day and final_day must be ordered, distinct integers")
    cohorts = config.get("cohorts", [])
    if not cohorts:
        raise ValueError("Bud input manifest has no cohorts")
    names = [entry["cohort"] for entry in cohorts]
    if len(names) != len(set(names)):
        raise ValueError("Each cohort name must be unique")
    measurements = []
    source_files = []
    for cohort in cohorts:
        sources = cohort.get("sources", [])
        if not sources:
            raise ValueError(f"Cohort {cohort['cohort']} has no sources")
        for source in sources:
            record_path = (config_path.parent / source["record"]).resolve()
            count_path = (config_path.parent / source["counts"]).resolve()
            decoded = _decode_pair(record_path, count_path)
            for role, path in (("filename_record", record_path), ("bud_counts", count_path)):
                source_files.append({"cohort": cohort["cohort"], "role": role, "path": str(path), "sha256": _sha256(path)})
            for row in decoded.to_dict("records"):
                parsed = _parse_name(row["original_file_name"], cohort, source)
                reason = ""
                if parsed["placeholder"]:
                    if not cohort.get("exclude_filename_placeholders", False):
                        raise ValueError(f"Found filename-marked placeholder {row['original_file_name']}; enable explicit placeholder exclusion")
                    reason = "filename-marked space holder"
                elif parsed["day"] not in (baseline_day, final_day):
                    reason = "outside selected baseline/final interval"
                measurements.append({
                    **row, **parsed, "cohort": cohort["cohort"],
                    "record_file": str(record_path), "count_file": str(count_path),
                    "included_in_endpoint": not bool(reason), "exclusion_reason": reason,
                })
    frame = pd.DataFrame(measurements)
    if frame.duplicated(["gland_id", "day"]).any():
        duplicates = frame.loc[frame.duplicated(["gland_id", "day"], keep=False), ["gland_id", "day", "original_file_name"]]
        raise ValueError(f"Duplicate measurements for the same gland/day:\n{duplicates.to_string(index=False)}")
    endpoints = []
    for gland_id, records in frame.groupby("gland_id", sort=True):
        selected = records[records.day.isin([baseline_day, final_day])]
        if selected.empty:
            continue
        if set(selected.day) != {baseline_day, final_day}:
            raise ValueError(f"{gland_id}: missing baseline or final observation")
        if selected.placeholder.any():
            if not selected.placeholder.all():
                raise ValueError(f"{gland_id}: placeholder status differs between endpoint days")
            continue
        if selected.dose_um.nunique() != 1:
            raise ValueError(f"{gland_id}: dose changes between endpoint days")
        baseline = selected.loc[selected.day.eq(baseline_day)].iloc[0]
        final = selected.loc[selected.day.eq(final_day)].iloc[0]
        if baseline.bud_number <= 0 or final.bud_number <= 0:
            raise ValueError(f"{gland_id}: log2 ratios require positive counts on both days; no automatic exclusion performed")
        ratio = int(final.bud_number) / int(baseline.bud_number)
        endpoints.append({
            "dose_um": float(baseline.dose_um), "value": math.log2(ratio),
            "cohort": baseline.cohort, "gland_id": gland_id,
            "group_id": baseline.group_id, "gland_number": baseline.gland_number,
            "baseline_day": baseline_day, "final_day": final_day,
            "baseline_buds": int(baseline.bud_number), "final_buds": int(final.bud_number),
            "ratio": ratio,
            "baseline_file": baseline.original_file_name, "final_file": final.original_file_name,
        })
    if not endpoints:
        raise ValueError("No complete, non-placeholder gland endpoints")
    result = pd.DataFrame(endpoints).sort_values(["dose_um", "cohort", "gland_id"]).reset_index(drop=True)
    audit = {
        "config_path": str(config_path), "config_sha256": _sha256(config_path),
        "schema_version": 1, "endpoint": f"log2(day{final_day} buds / day{baseline_day} buds)",
        "pairing": "Explicit cohort, group and gland number; filename-record joins are one-to-one.",
        "exclusion_policy": "Only explicitly enabled filename-marked placeholders and acquisition days outside the selected interval; no count threshold filter.",
        "input_measurements": len(frame), "gland_endpoints": len(result),
        "source_files": source_files,
        "measurements": frame.to_dict("records"),
        "excluded_measurements": frame.loc[~frame.included_in_endpoint].to_dict("records"),
        "n_by_dose": [{"dose_um": float(dose), "n": int(n)} for dose, n in result.groupby("dose_um").size().items()],
        "n_by_cohort_and_dose": result.groupby(["cohort", "dose_um"]).size().rename("n").reset_index().to_dict("records"),
    }
    return result, audit
