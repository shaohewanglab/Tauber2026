import unittest

import numpy as np
from scipy import stats

from trend_statistics import permutation_spearman


class PermutationSpearmanTests(unittest.TestCase):
    def test_exact_tied_case_matches_scipy_independent_reference(self):
        x = np.array([1, 1, 2, 2, 3, 3])
        y = np.array([6, 6, 4, 3, 1, 1])
        result = permutation_spearman(x, y)
        reference = stats.permutation_test(
            (y,), lambda v: abs(stats.spearmanr(x, v).statistic),
            permutation_type="pairings", alternative="greater", n_resamples=np.inf)
        self.assertAlmostEqual(result["rho"], stats.spearmanr(x, y).statistic)
        self.assertAlmostEqual(result["p_raw"], reference.pvalue)
        self.assertEqual(result["permutation_mode"], "exhaustive")

    def test_monte_carlo_is_reproducible_and_does_not_return_zero(self):
        x = np.arange(12)
        first = permutation_spearman(x, -x, n_resamples=999, seed=7, batch=100)
        second = permutation_spearman(x, -x, n_resamples=999, seed=7, batch=100)
        self.assertEqual(first, second)
        self.assertAlmostEqual(first["rho"], -1)
        self.assertEqual(first["p_raw"], .001)
        self.assertEqual(first["extreme_permutations"], 0)

    def test_rank_statistic_unchanged_by_monotonic_transform(self):
        x = np.arange(6)
        y = [1, 2, 2, 5, 4, 7]
        first = permutation_spearman(x, y)
        second = permutation_spearman(x, np.log2(y))
        self.assertEqual(first, second)

    def test_constant_or_missing_data_are_not_tested(self):
        for x, y in [([1, 1, 1], [2, 3, 4]), ([1, 2, 3], [2, 2, 2]),
                     ([1, 2, 3], [2, np.nan, 4])]:
            with self.assertRaises(ValueError):
                permutation_spearman(x, y)


if __name__ == "__main__":
    unittest.main()
