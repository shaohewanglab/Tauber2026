"""Tests for explicit gland pairing and rejection of ambiguous count inputs."""
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from titration_bud_io import load_bud_data


class BudInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.config = {
            "schema_version": 1,
            "cohorts": [{
                "cohort": "test", "filename_format": "date_group_gland",
                "date_to_day": {"250101": 0, "250102": 1},
                "dose_um_by_group": {"1": 0},
                "sources": [{"record": "keys.tsv", "counts": "counts.tsv"}],
            }],
        }
        self.keys = pd.DataFrame({
            "original_file_name": ["250101-G01_0000.TIFF", "250102-G01_0000.TIFF"],
            "converted_file_name": ["a.TIFF", "b.TIFF"],
        })
        self.counts = pd.DataFrame({"file_name": ["b.TIFF", "a.TIFF"], "bud_number": [8, 2]})

    def load(self):
        config_path = self.folder / "config.json"
        config_path.write_text(json.dumps(self.config))
        self.keys.to_csv(self.folder / "keys.tsv", sep="\t", index=False)
        self.counts.to_csv(self.folder / "counts.tsv", sep="\t", index=False)
        return load_bud_data(config_path)

    def test_pairing_ignores_row_order(self):
        self.assertEqual(self.load()[0].value.iloc[0], 2)

    def test_one_bud_is_not_excluded(self):
        self.counts.bud_number = [8, 1]
        self.assertEqual(self.load()[0].value.iloc[0], 3)

    def test_missing_count_is_rejected(self):
        self.counts = self.counts.iloc[:1]
        with self.assertRaisesRegex(ValueError, "mismatch"):
            self.load()

    def test_duplicate_count_is_rejected(self):
        self.counts = pd.concat([self.counts, self.counts.iloc[:1]])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.load()

    def test_zero_count_is_rejected(self):
        self.counts.bud_number = [0, 2]
        with self.assertRaisesRegex(ValueError, "positive"):
            self.load()

    def test_fractional_count_is_rejected(self):
        self.counts.bud_number = [1.5, 2]
        with self.assertRaisesRegex(ValueError, "integers"):
            self.load()

    def test_incomplete_gland_pair_is_rejected(self):
        self.keys.original_file_name = ["250101-G01_0000.TIFF", "250102-G01_0001.TIFF"]
        with self.assertRaisesRegex(ValueError, "missing baseline"):
            self.load()

    def test_unknown_date_is_rejected(self):
        self.keys.original_file_name = ["250101-G01_0000.TIFF", "250103-G01_0000.TIFF"]
        with self.assertRaisesRegex(ValueError, "absent from date"):
            self.load()

    def test_inconsistent_placeholder_is_rejected(self):
        self.config["cohorts"][0]["exclude_filename_placeholders"] = True
        self.keys.original_file_name = ["250101-G01_0000-space-holder.TIFF", "250102-G01_0000.TIFF"]
        with self.assertRaisesRegex(ValueError, "placeholder status differs"):
            self.load()


if __name__ == "__main__":
    unittest.main()
