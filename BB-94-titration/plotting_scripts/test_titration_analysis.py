"""Checks on inferential direction, variance handling and multiplicity."""
import unittest

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from titration_analysis import analyze, significance


def endpoint_data(groups):
    return pd.DataFrame([
        {"dose_um": dose, "value": float(value), "cohort": "a", "gland_id": f"{dose}-{index}"}
        for dose, values in enumerate(groups) for index, value in enumerate(values)
    ])


class StatisticsTests(unittest.TestCase):
    def setUp(self):
        self.groups = [np.array([1., 1.2, .8, 1.4, .9]), np.array([2., 4., 3., 7., 5., 2.5])]
        self.data = pd.DataFrame([{"dose_um": dose, "value": value, "cohort": "a", "gland_id": f"{dose}-{i}"}
                                  for dose, values in enumerate(self.groups) for i, value in enumerate(values)])

    def test_default_is_dunnett_vs_vehicle_with_pooled_t_for_two_groups(self):
        _, actual, diagnostics = analyze(self.data)
        reference = stats.ttest_ind(self.groups[1], self.groups[0], equal_var=True)
        self.assertEqual(len(actual), 1)
        self.assertIn("Dunnett", actual.method.iloc[0])
        self.assertAlmostEqual(actual.t_b.iloc[0], reference.statistic)
        self.assertAlmostEqual(actual.p_adjusted.iloc[0], reference.pvalue)
        self.assertEqual(diagnostics["df"], len(self.data)-2)

    def test_two_group_games_howell_equals_welch(self):
        _, result, _ = analyze(self.data, "games-howell")
        reference = stats.ttest_ind(self.groups[1], self.groups[0], equal_var=False)
        ci = reference.confidence_interval()
        self.assertAlmostEqual(result.p_adjusted.iloc[0], reference.pvalue, places=11)
        self.assertAlmostEqual(result.ci95_low.iloc[0], ci.low, places=9)
        self.assertAlmostEqual(result.ci95_high.iloc[0], ci.high, places=9)
        self.assertGreater(result.difference_b_minus_a.iloc[0], 0)

    def test_tukey_matches_scipy_and_controls_family(self):
        extra = pd.DataFrame({"dose_um": 2, "value": [4., 2., 3., 5., 6.], "cohort": "a", "gland_id": [f"x{i}" for i in range(5)]})
        data = pd.concat([self.data, extra], ignore_index=True)
        _, actual, _ = analyze(data, "tukey")
        reference = stats.tukey_hsd(*self.groups, extra.value.to_numpy())
        self.assertEqual(len(actual), 3)
        for row in actual.itertuples():
            i, j = int(row.dose_a_um), int(row.dose_b_um)
            self.assertAlmostEqual(row.p_adjusted, reference.pvalue[j, i], places=13)
            self.assertAlmostEqual(row.difference_b_minus_a, reference.statistic[j, i], places=13)

    def test_normalization_cannot_change_significance(self):
        shifted = self.data.copy()
        shifted["value"] = shifted.value*17 + 9
        _, before, _ = analyze(self.data, "games-howell")
        _, after, _ = analyze(shifted, "games-howell")
        np.testing.assert_allclose(before.p_adjusted, after.p_adjusted, atol=1e-13)
        np.testing.assert_allclose(before.difference_b_minus_a*17, after.difference_b_minus_a)

    def test_duplicate_glands_rejected(self):
        with self.assertRaisesRegex(ValueError, "exactly one"):
            analyze(pd.concat([self.data, self.data.iloc[:1]]))

    def test_star_thresholds(self):
        self.assertEqual([significance(p) for p in [1, .05, .049, .009, .0009, .00009]], ["ns", "ns", "*", "**", "***", "****"])
        with self.assertRaises(ValueError):
            significance(np.nan)


class MannWhitneyTests(unittest.TestCase):
    def test_explicit_mann_whitney_exact_seven_per_group_and_effect_direction(self):
        data = endpoint_data([np.arange(7), np.arange(7, 14)])
        _, result, _ = analyze(data, "mann-whitney")
        row = result.iloc[0]
        # There are C(14, 7)=3432 assignments; the two completely separated
        # assignments form the two-sided exact tail.
        self.assertAlmostEqual(row.p_raw, 2 / 3432, places=14)
        self.assertAlmostEqual(row.p_adjusted, 2 / 3432, places=14)
        self.assertEqual(row.U_b, 49)
        self.assertEqual(row.rank_biserial_b_minus_a, 1)
        self.assertEqual(row.mean_difference_b_minus_a, 7)
        self.assertIn("exact", row.p_method.lower())
        self.assertEqual(row.annotation, "***")
        self.assertNotIn("ci95_low", result.columns)
        self.assertNotIn("ci95_high", result.columns)

    def test_ties_use_asymptotic_tie_and_continuity_correction(self):
        a = np.array([0., 0., 1., 1., 2., 2., 3.])
        b = np.array([1., 1., 2., 3., 3., 4., 4.])
        _, result, _ = analyze(endpoint_data([a, b]), "mann-whitney")
        reference = stats.mannwhitneyu(b, a, alternative="two-sided", method="asymptotic", use_continuity=True)
        row = result.iloc[0]
        self.assertAlmostEqual(row.U_b, reference.statistic, places=14)
        self.assertAlmostEqual(row.p_raw, reference.pvalue, places=14)
        self.assertIn("asymptotic", row.p_method.lower())
        self.assertGreater(row.rank_biserial_b_minus_a, 0)

    def test_nine_per_group_without_ties_uses_asymptotic(self):
        a, b = np.arange(9), np.arange(9, 18)
        _, result, _ = analyze(endpoint_data([a, b]), "mann-whitney")
        reference = stats.mannwhitneyu(b, a, alternative="two-sided", method="auto", use_continuity=True)
        self.assertAlmostEqual(result.p_raw.iloc[0], reference.pvalue, places=14)
        self.assertIn("asymptotic", result.p_method.iloc[0].lower())

    def test_holm_family_includes_six_control_comparisons_only(self):
        groups = [np.arange(7) + 10 * dose for dose in range(7)]
        _, result, _ = analyze(endpoint_data(groups), "mann-whitney")
        self.assertEqual(len(result), 6)
        self.assertTrue(result.dose_a_um.eq(0).all())
        self.assertEqual(set(result.dose_b_um), set(range(1, 7)))
        expected = multipletests(result.p_raw.to_numpy(), method="holm")[1]
        np.testing.assert_allclose(result.p_adjusted, expected, rtol=1e-13)
        np.testing.assert_allclose(result.p_adjusted, 6 * 2 / 3432, rtol=1e-13)
        self.assertTrue(result.annotation.eq("**").all())

    def test_treatment_pair_separation_does_not_expand_correction_family(self):
        # Both alternatives completely separate every treatment from control.
        # Moving the treatment groups apart changes treatment-versus-treatment
        # comparisons, which must not enter the requested control-only family.
        overlapping = [np.arange(7)] + [np.arange(10, 17)] * 6
        separated = [np.arange(7) + 10 * dose for dose in range(7)]
        _, before, _ = analyze(endpoint_data(overlapping), "mann-whitney")
        _, after, _ = analyze(endpoint_data(separated), "mann-whitney")
        np.testing.assert_allclose(before.p_raw, after.p_raw, rtol=1e-13)
        np.testing.assert_allclose(before.p_adjusted, after.p_adjusted, rtol=1e-13)

    def test_all_identical_observations_are_valid_and_not_significant(self):
        _, result, _ = analyze(endpoint_data([np.ones(7), np.ones(7)]), "mann-whitney")
        row = result.iloc[0]
        self.assertEqual(row.p_raw, 1)
        self.assertEqual(row.p_adjusted, 1)
        self.assertEqual(row.U_b, 24.5)
        self.assertEqual(row.rank_biserial_b_minus_a, 0)
        self.assertEqual(row.annotation, "ns")

    def test_rank_result_is_unchanged_by_strictly_increasing_transform(self):
        data = endpoint_data([np.arange(1, 8), np.arange(5, 12)])
        transformed = data.copy()
        transformed["value"] = np.log2(transformed.value)
        _, original, _ = analyze(data, "mann-whitney")
        _, changed, _ = analyze(transformed, "mann-whitney")
        np.testing.assert_allclose(original.p_adjusted, changed.p_adjusted, rtol=1e-13)
        np.testing.assert_allclose(original.rank_biserial_b_minus_a, changed.rank_biserial_b_minus_a, rtol=1e-13)


if __name__ == "__main__":
    unittest.main()
