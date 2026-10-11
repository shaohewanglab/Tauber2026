#!/usr/bin/env python3
"""Check source integrity and reproduce the accepted numerical results."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def main():
    sources = json.loads((ROOT/"source-data-manifest.json").read_text())
    for source in sources:
        path = ROOT/source["path"]
        assert path.is_file(), f"Missing source: {source['path']}"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == source["sha256"], f"Source changed: {source['path']}"
    reference = json.loads((ROOT/"validation/expected-results.json").read_text())
    for specification in reference["tables"]:
        path = ROOT/specification["path"]
        expected = pd.DataFrame(specification["records"])
        observed = pd.read_csv(path, dtype=specification.get("dtype", {}))
        columns = specification["columns"]
        sort_by = specification["sort_by"]
        expected = expected[columns].sort_values(sort_by).reset_index(drop=True)
        observed = observed[columns].sort_values(sort_by).reset_index(drop=True)
        pd.testing.assert_frame_equal(observed, expected, check_dtype=False, check_exact=False, rtol=1e-10, atol=1e-12, obj=specification["path"])
    for plot in reference["plots"]:
        for extension in ["pdf", "svg", "png"]:
            path = ROOT/(plot+"."+extension)
            assert path.is_file() and path.stat().st_size > 0, f"Missing plot: {path}"
    report = {"source_tables_verified": len(sources), "numerical_tables_verified": len(reference["tables"]),
              "plots_verified": len(reference["plots"]), "relative_tolerance": 1e-10, "absolute_tolerance": 1e-12,
              "status": "passed", "note": "Expected results were captured from the accepted manuscript analyses; statistical unit, numerical endpoints, SD and adjusted p-values are checked. Font substitution may change appearance only."}
    (ROOT/"validation/last-verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(f"Verified {len(sources)} source tables, {len(reference['tables'])} numerical tables and {len(reference['plots'])} plots.")


if __name__ == "__main__":
    main()
