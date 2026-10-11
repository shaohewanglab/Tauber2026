"""Tests for biological replicate aggregation and malformed-input rejection."""
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from titration_collagen_io import load_collagen_data


class CollagenInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.control = pd.DataFrame({
            "file_name": ["1", "1", "2"], "id": ["a", "b", "a"],
            "channel": [2, 2, 2], "AUC": [1, 3, 6],
        })
        self.treated = pd.DataFrame({
            "file_name": ["1", "1"], "id": ["a", "b"],
            "channel": [2, 2], "AUC": [6, 10],
        })
        self.config = {
            "schema_version": 1,
            "inputs": [
                {"path": "control.csv", "cohort": "test", "dose_um": 0},
                {"path": "treated.csv", "cohort": "test", "dose_um": 5},
            ],
        }

    def load(self):
        self.control.to_csv(self.folder / "control.csv", index=False)
        self.treated.to_csv(self.folder / "treated.csv", index=False)
        config_path = self.folder / "config.json"
        config_path.write_text(json.dumps(self.config))
        return load_collagen_data(config_path)

    def test_equal_gland_weighting_with_unequal_scan_counts(self):
        data, audit = self.load()
        self.assertEqual(data["raw_auc"].tolist(), [2, 6, 8])
        self.assertEqual(data["value"].tolist(), [0.5, 1.5, 2])
        self.assertEqual(data["n_scans"].tolist(), [2, 1, 2])
        self.assertEqual(audit["normalization"]["control_gland_mean_auc"], 4)
        self.assertAlmostEqual(audit["normalization"]["notebook_pooled_control_scan_mean_auc"], 10 / 3)
        self.assertEqual(data["gland_id"].nunique(), 3)

    def test_duplicate_scan_is_rejected(self):
        self.control.loc[1, "id"] = "a"
        with self.assertRaisesRegex(ValueError, "Duplicate gland/scan"):
            self.load()

    def test_nonfinite_auc_is_rejected(self):
        self.control["AUC"] = self.control["AUC"].astype(float)
        self.control.loc[1, "AUC"] = float("inf")
        with self.assertRaisesRegex(ValueError, "finite"):
            self.load()

    def test_channel_mix_is_rejected(self):
        self.control.loc[1, "channel"] = 3
        with self.assertRaisesRegex(ValueError, "channel 2"):
            self.load()

    def test_missing_gland_id_is_rejected(self):
        self.control.loc[1, "file_name"] = " "
        with self.assertRaisesRegex(ValueError, "gland IDs"):
            self.load()

    def test_too_few_scans_is_rejected(self):
        self.config["minimum_scans_per_gland"] = 2
        with self.assertRaisesRegex(ValueError, "fewer than 2"):
            self.load()

    def test_wrong_expected_gland_count_is_rejected(self):
        self.config["inputs"][0]["expected_glands"] = 3
        with self.assertRaisesRegex(ValueError, "expected 3 glands"):
            self.load()

    def test_duplicate_source_is_rejected(self):
        self.config["inputs"].append(self.config["inputs"][0].copy())
        with self.assertRaisesRegex(ValueError, "listed twice"):
            self.load()


if __name__ == "__main__":
    unittest.main()
