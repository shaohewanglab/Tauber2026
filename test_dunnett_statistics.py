import unittest

import numpy as np
from scipy import stats

from dunnett_statistics import dunnett_vs_control, _adjusted_tail


class DunnettTests(unittest.TestCase):
    def test_one_comparison_is_pooled_two_sample_t(self):
        control, sample = [1, 2, 4], [2, 4, 5, 7, 8]
        actual = dunnett_vs_control([sample], control)
        reference = stats.ttest_ind(sample, control, equal_var=True)
        np.testing.assert_allclose(actual['statistic'], [reference.statistic], rtol=1e-13)
        np.testing.assert_allclose(actual['p_adjusted'], [reference.pvalue], rtol=1e-13)

    def test_unequal_sizes_match_scipy_documented_example(self):
        result = dunnett_vs_control([[3.8, 2.7, 4, 2.4], [2.8, 3.4, 3.7, 2.2, 2]],
                                    [2.9, 3, 2.5, 2.6, 3.2])
        np.testing.assert_allclose(result['statistic'], [.90874545, -.05007117], atol=1e-8)
        # Published SciPy values are from randomized multivariate-t integration.
        np.testing.assert_allclose(result['p_adjusted'], [.58325114, .99819341], atol=2e-4)

    def test_small_model_tails_converge_and_identical_contrasts_match(self):
        samples = [[5, 6, 5], [4, 3, 4], [1, 1, 3], [1, 1, 1], [1, 1, 1]]
        result = dunnett_vs_control(samples, [6, 6, 6])
        coarse = [_adjusted_tail(float(abs(t)), (3,)*5, 3, 12, nodes=128)
                  for t in result['statistic']]
        np.testing.assert_allclose(result['p_adjusted'], coarse, rtol=1e-8, atol=1e-16)
        self.assertEqual(result['p_adjusted'][-1], result['p_adjusted'][-2])
        self.assertEqual(result['df'], 12)
        self.assertAlmostEqual(result['pooled_variance'], 1/3)
        self.assertTrue(np.all(result['p_adjusted'] >= result['p_raw']))
        self.assertTrue(np.all(result['p_adjusted'] <= 5*result['p_raw']))

    def test_affine_scale_and_group_order_do_not_change_p(self):
        control = np.array([1., 2, 4, 6])
        groups = [np.array([2., 3, 6]), np.array([4., 5, 6, 8, 9])]
        a = dunnett_vs_control(groups, control)
        b = dunnett_vs_control([g*3+7 for g in reversed(groups)], control*3+7)
        np.testing.assert_allclose(a['p_adjusted'], b['p_adjusted'][::-1], rtol=1e-10)

    def test_invalid_and_degenerate_inputs_rejected(self):
        for samples, control in [([[2]], [1, 2]), ([[2, np.nan]], [1, 2]),
                                 ([[1, 1]], [2, 2]), ([], [1, 2])]:
            with self.assertRaises(ValueError):
                dunnett_vs_control(samples, control)


if __name__ == '__main__':
    unittest.main()
