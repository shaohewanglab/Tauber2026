"""Tests that check cohort adjustment, not merely code execution."""

import unittest

import numpy as np
import pandas as pd

from titration_sensitivity import bud_cohort_sensitivity


class CohortSensitivityTests(unittest.TestCase):
    def test_contrasts_are_invariant_to_additive_cohort_shifts_and_names(self):
        # Every cohort contains the same centered noise and known dose effects.
        rows = []
        for cohort in ["A", "B", "C"]:
            for dose, effect in [(0.0, 0.0), (0.5, -0.75), (5.0, -1.5)]:
                for gland, noise in enumerate([-0.3, -0.1, 0.1, 0.3]):
                    rows.append({
                        "cohort": cohort, "dose_um": dose,
                        "gland_id": f"{cohort}:{dose}:{gland}",
                        "value": 2.0 + effect + noise,
                    })
        original = pd.DataFrame(rows)
        shifted = original.copy()
        shifted["value"] += shifted["cohort"].map({"A": 10.0, "B": -7.0, "C": 3.5})
        shifted["cohort"] = shifted["cohort"].map({"A": "z", "B": "first", "C": "middle"})
        shifted = shifted.sample(frac=1, random_state=13)
        first, metadata = bud_cohort_sensitivity(original)
        second, _ = bud_cohort_sensitivity(shifted)
        np.testing.assert_allclose(first["mean_difference"], [-0.75, -1.5], atol=1e-12)
        for column in ["mean_difference", "standard_error", "p_unadjusted", "p_adjusted", "ci95_low", "ci95_high"]:
            np.testing.assert_allclose(first[column], second[column], rtol=1e-10, atol=1e-12)
        self.assertTrue(metadata["use_t"])
        self.assertEqual(metadata["n_comparisons"], 2)
        self.assertEqual(metadata["residual_df"], 31)

    def test_repeated_gland_is_rejected(self):
        data = pd.DataFrame({"dose_um": [0, 1], "value": [1, 2], "cohort": ["A", "B"], "gland_id": ["same", "same"]})
        with self.assertRaisesRegex(ValueError, "one response per gland"):
            bud_cohort_sensitivity(data)


if __name__ == "__main__":
    unittest.main()
