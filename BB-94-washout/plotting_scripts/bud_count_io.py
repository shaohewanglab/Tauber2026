"""Strict decoding and metadata validation for randomized bud-count files.

``load_pairs([(record_path, count_path), ...])`` returns one row per image.
Filename dates are acquisition dates, not gland identifiers: longitudinal
subjects are identified by experiment, treatment group, and gland number.
Only NumPy and pandas are required beyond the Python standard library.
"""

import csv
from decimal import Decimal, InvalidOperation
from pathlib import Path
import re
from urllib.parse import quote

import numpy as np
import pandas as pd


OUTPUT_COLUMNS = [
    "experiment", "group", "gland", "time_h", "bud_number",
    "original_file_name", "converted_file_name", "record_file", "count_file",
    "subject_id", "plate", "time_unit", "time_value",
]
_FILENAME = re.compile(
    r"^(?:\d{6}|\d{8})-(?:SMG|plate(?P<plate>\d+))-"
    r"(?P<group>.+?)-(?P<unit>day|d|t)(?P<time>\d+(?:\.\d+)?)_"
    r"(?P<gland>[^.]+)\.[^.]+$",
    re.IGNORECASE,
)


class InputError(ValueError):
    """An input cannot be decoded safely or is internally inconsistent."""


def _short(values):
    values = sorted(str(value) for value in values)
    return ", ".join(repr(value) for value in values[:5]) + (
        " ..." if len(values) > 5 else ""
    )


def _read_table(path, required):
    path = Path(path).expanduser().resolve()
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            sample = stream.read(8192)
            stream.seek(0)
            if not sample.strip():
                raise InputError("Input is empty: {}".format(path))
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
            except csv.Error as exc:
                raise InputError(
                    "Could not identify tab or comma delimiters in {}".format(path)
                ) from exc
            reader = csv.reader(stream, dialect)
            rows = [row for row in reader if any(cell.strip() for cell in row)]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise InputError("Could not read {}: {}".format(path, exc)) from exc

    headers = [cell.strip() for cell in rows[0]]
    if len(set(headers)) != len(headers) or any(not name for name in headers):
        raise InputError("Blank or duplicate column headers in {}".format(path))
    missing = set(required) - set(headers)
    if missing:
        raise InputError("{} is missing columns: {}".format(path, _short(missing)))
    if len(rows) == 1:
        raise InputError("Input has a header but no data: {}".format(path))
    for index, row in enumerate(rows[1:], 2):
        if len(row) != len(headers):
            raise InputError(
                "{}: data row {} has {} fields; expected {}".format(
                    path, index, len(row), len(headers)
                )
            )
    return pd.DataFrame(rows[1:], columns=headers), path


def _names(series, path, column):
    """Normalize surrounding whitespace and POSIX/Windows path prefixes."""
    normalized = series.map(lambda value: value.strip().replace("\\", "/").rsplit("/", 1)[-1])
    if normalized.eq("").any():
        raise InputError("{} has blank {} values".format(path, column))
    duplicates = normalized[normalized.duplicated(keep=False)]
    if not duplicates.empty:
        raise InputError(
            "{} has duplicate {} values after basename/whitespace normalization: {}".format(
                path, column, _short(set(duplicates))
            )
        )
    return normalized


def _counts(series, path):
    result = []
    for index, text in enumerate(series, 2):
        try:
            number = Decimal(text.strip())
            valid = (
                number.is_finite() and number >= 0
                and number == number.to_integral_value()
                and number <= np.iinfo(np.int64).max
            )
        except InvalidOperation:
            valid = False
        if not valid:
            raise InputError(
                "{}: bud_number on row {} must be a finite nonnegative integer; got {!r}".format(
                    path, index, text
                )
            )
        result.append(int(number))
    return np.asarray(result, dtype=np.int64)


def _labels(series, label, path):
    values = series.map(lambda value: str(value).strip())
    if values.eq("").any():
        raise InputError("{} has blank {} values".format(path, label))
    return values


def _natural(value):
    return tuple(
        (1, int(piece)) if piece.isdigit() else (0, piece.casefold())
        for piece in re.split(r"(\d+)", str(value))
    )


def load_pairs(pairs, experiment="experiment1", metadata_path=None):
    """Decode record/count pairs into a validated, deterministically sorted table.

    Record tables require ``original_file_name, converted_file_name``; count
    tables require ``file_name, bud_number``. Comma and tab delimiters are
    detected automatically. Each pair must match exactly, with no duplicate,
    missing, or extra image names. Multiple pairs may cover different days of
    one experiment; they must not repeat an image or a gland/time observation.

    By default, names such as ``210325-SMG-ctrl-t00_1.nd2`` or
    ``210401-plate2-control-day1_5.nd2`` supply treatment, time, and gland.
    Six- and eight-digit dates, decimal times, and either SMG or plateN are
    supported; d0/d1/d2 are aliases of day0/day1/day2. Structural tokens are
    case-insensitive. The known group aliases
    ctrl/control and BB94/BB-94 become ctrl and BB94, respectively. Day values
    are converted to hours in time_h, while time_unit and time_value retain
    their source units. Plates are included in gland identity (plate2_5), so
    glands numbered identically on different plates remain distinct. SMG-only
    filenames retain their original gland identifiers.

    A metadata table instead supplies those fields
    through ``original_file_name, group, gland, time_h`` plus an optional
    ``experiment`` column. It must map every decoded image exactly once and
    contain no unused rows. Experiment labels distinguish independent cohorts.
    Metadata is authoritative: its group/gland labels are preserved, plate is
    empty, and time_unit/time_value are hour/time_h. Additional unit or plate
    columns in metadata are not used; include plate identity in gland yourself.

    Output columns, in order, are ``OUTPUT_COLUMNS``. ``subject_id`` joins
    URL-escaped experiment, group, and gland with ``|``. Acquisition date is
    deliberately excluded so a gland keeps its identity across time points.
    Raises InputError (a ValueError subclass) for invalid inputs.
    """
    pairs = list(pairs)
    if not pairs:
        raise InputError("At least one record/count pair is required")
    experiment = str(experiment).strip()
    if not experiment:
        raise InputError("Experiment label must not be blank")
    decoded = []
    for pair in pairs:
        if len(pair) != 2:
            raise InputError("Each pair must contain a record path and a count path")
        record, record_path = _read_table(
            pair[0], ["original_file_name", "converted_file_name"]
        )
        counts, count_path = _read_table(pair[1], ["file_name", "bud_number"])
        for name in ("original_file_name", "converted_file_name"):
            record[name] = _names(record[name], record_path, name)
        counts["file_name"] = _names(counts["file_name"], count_path, "file_name")
        expected = set(record["converted_file_name"])
        observed = set(counts["file_name"])
        if expected != observed:
            details = []
            if expected - observed:
                details.append("missing counts for " + _short(expected - observed))
            if observed - expected:
                details.append("counts without a record for " + _short(observed - expected))
            raise InputError(
                "Record/count pair does not match ({}; {}): {}".format(
                    record_path, count_path, "; ".join(details)
                )
            )
        counts["bud_number"] = _counts(counts["bud_number"], count_path)
        frame = record[["original_file_name", "converted_file_name"]].merge(
            counts[["file_name", "bud_number"]],
            left_on="converted_file_name", right_on="file_name", validate="one_to_one"
        ).drop(columns="file_name")
        frame["record_file"] = str(record_path)
        frame["count_file"] = str(count_path)
        decoded.append(frame)

    result = pd.concat(decoded, ignore_index=True)
    for column in ("original_file_name", "converted_file_name"):
        duplicates = result.loc[result[column].duplicated(keep=False), column]
        if not duplicates.empty:
            raise InputError("Repeated {} across pairs: {}".format(column, _short(set(duplicates))))

    if metadata_path is not None:
        metadata, metadata_path = _read_table(
            metadata_path, ["original_file_name", "group", "gland", "time_h"]
        )
        metadata["original_file_name"] = _names(
            metadata["original_file_name"], metadata_path, "original_file_name"
        )
        expected = set(result["original_file_name"])
        supplied = set(metadata["original_file_name"])
        if expected != supplied:
            details = []
            if expected - supplied:
                details.append("missing mappings for " + _short(expected - supplied))
            if supplied - expected:
                details.append("unused mappings for " + _short(supplied - expected))
            raise InputError("Metadata coverage does not match inputs: " + "; ".join(details))
        if "experiment" not in metadata:
            metadata["experiment"] = experiment
        for label in ("experiment", "group", "gland"):
            metadata[label] = _labels(metadata[label], label, metadata_path)
        metadata["time_h"] = pd.to_numeric(metadata["time_h"].str.strip(), errors="coerce")
        if not np.isfinite(metadata["time_h"].to_numpy(dtype=float)).all():
            raise InputError("{}: time_h must contain finite numbers".format(metadata_path))
        result = result.merge(
            metadata[["original_file_name", "experiment", "group", "gland", "time_h"]],
            on="original_file_name", validate="one_to_one"
        )
        result["plate"] = ""
        result["time_unit"] = "hour"
        result["time_value"] = result["time_h"]
    else:
        parsed = []
        for name in result["original_file_name"]:
            match = _FILENAME.fullmatch(name)
            if match is None:
                raise InputError(
                    "Cannot parse {!r}; expected date-SMG-group-tHH_gland.ext "
                    "or date-plateN-group-dayD_gland.ext, with a six- or "
                    "eight-digit date and a nonnegative time "
                    "(for example 210325-SMG-ctrl-t00_1.nd2 or "
                    "210401-plate2-control-day1_5.nd2). "
                    "Use --metadata for other filename conventions.".format(name)
                )
            group = match["group"].strip()
            group = {"ctrl": "ctrl", "control": "ctrl", "bb94": "BB94", "bb-94": "BB94"}.get(group.casefold(), group)
            plate = "plate{}".format(int(match["plate"])) if match["plate"] is not None else ""
            gland = "{}_{}".format(plate, match["gland"]) if plate else match["gland"]
            unit = "day" if match["unit"].casefold() in ("day", "d") else "hour"
            time_value = float(match["time"])
            time_h = time_value * 24 if unit == "day" else time_value
            parsed.append((experiment, group, gland, time_h, plate, unit, time_value))
        result[["experiment", "group", "gland", "time_h", "plate", "time_unit", "time_value"]] = (
            pd.DataFrame(parsed, index=result.index)
        )

    for label in ("experiment", "group", "gland"):
        result[label] = _labels(result[label], label, "Decoded metadata")
    if not np.isfinite(result["time_h"].to_numpy(dtype=float)).all():
        raise InputError("Decoded time_h must contain finite numbers")

    cells = ["experiment", "group", "gland", "time_h"]
    duplicates = result.loc[result.duplicated(cells, keep=False), cells]
    if not duplicates.empty:
        raise InputError(
            "Duplicate experiment/group/gland/time observations: "
            + _short(tuple(row) for row in duplicates.itertuples(index=False, name=None))
        )
    result["subject_id"] = [
        "|".join(quote(str(value), safe="") for value in row)
        for row in result[["experiment", "group", "gland"]].itertuples(index=False, name=None)
    ]
    order = sorted(
        result.index,
        key=lambda index: (
            _natural(result.at[index, "experiment"]),
            _natural(result.at[index, "group"]),
            _natural(result.at[index, "gland"]),
            float(result.at[index, "time_h"]),
            result.at[index, "original_file_name"],
        )
    )
    return result.loc[order, OUTPUT_COLUMNS].reset_index(drop=True)
