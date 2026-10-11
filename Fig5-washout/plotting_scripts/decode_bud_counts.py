#!/usr/bin/env python3
"""Decode randomized bud counts without changing the source files.

Example::

    python decode_bud_counts.py --pair fileNameRecord.txt budCount.txt \
        --experiment washout-2021 --output decoded-bud-counts.csv
"""

import argparse
from pathlib import Path

from bud_count_io import InputError, load_pairs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--pair", nargs=2, action="append", required=True, metavar=("RECORD", "COUNTS"),
        help="Paired filename record and bud-count table; repeat for additional days."
    )
    parser.add_argument(
        "--experiment", default="experiment1",
        help="Cohort identifier (default: experiment1); dates do not identify glands."
    )
    parser.add_argument(
        "--metadata", metavar="CSV_OR_TSV",
        help="Optional original_file_name,group,gland,time_h table; experiment is optional."
    )
    parser.add_argument(
        "--output", type=Path, default=Path("decoded-bud-counts.csv"),
        help="Output CSV (default: decoded-bud-counts.csv in the current folder)."
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing output CSV; never overwrite an input.")
    args = parser.parse_args(argv)
    output = args.output.expanduser().resolve()
    inputs = [Path(path).expanduser().resolve() for pair in args.pair for path in pair]
    if args.metadata:
        inputs.append(Path(args.metadata).expanduser().resolve())
    if output in inputs or (output.exists() and any(output.samefile(path) for path in inputs if path.exists())):
        parser.error("Output must not overwrite any input file: {}".format(output))
    if output.exists() and not args.overwrite:
        parser.error("Output already exists: {}. Use --overwrite to replace it.".format(output))
    try:
        decoded = load_pairs(args.pair, experiment=args.experiment, metadata_path=args.metadata)
        output.parent.mkdir(parents=True, exist_ok=True)
        decoded.to_csv(output, index=False, mode="w" if args.overwrite else "x", lineterminator="\n")
    except (InputError, OSError) as exc:
        parser.error(str(exc))
    print("Decoded {} images from {} glands into {}".format(len(decoded), decoded["subject_id"].nunique(), output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
