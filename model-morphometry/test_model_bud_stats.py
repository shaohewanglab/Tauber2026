import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from model_bud_stats import analyze_bud_counts


def observations(groups):
    return pd.DataFrame([{"z": z, "replicate": i+1, "N_buds": value}
                         for z, values in groups.items() for i, value in enumerate(values)])


class ModelBudStatisticsTests(unittest.TestCase):
    def test_complete_shared_counts_use_five_dunnett_comparisons(self):
        source = pd.read_csv(Path(__file__).parent/"ModelData/bud_counts_summary.csv")
        data = observations({row.z: [row.Nbud1, row.Nbud2, row.Nbud3]
                             for row in source.itertuples()})
        result = analyze_bud_counts(data, source.z.tolist())
        self.assertEqual(len(result), 5)
        self.assertTrue((result.n_baseline == 3).all())
        self.assertTrue((result.n_comparison == 3).all())
        self.assertTrue((result.status == "tested").all())
        # Independently reviewed equal-n multivariate-t quadrature reference.
        np.testing.assert_allclose(result.p_adjusted,
            [.509792084207812, .00141820900225409, 3.91365176197465e-6,
             8.41973980730576e-7, 8.41973980730576e-7], rtol=1e-8)
        self.assertEqual(result.annotation.tolist(), ["ns", "**", "****", "****", "****"])
        self.assertTrue((result.df == 12).all())
        np.testing.assert_allclose(result.pooled_variance, 1/3)
        self.assertEqual(len(data), 18)

    def test_incomplete_family_is_unavailable_not_nonsignificant(self):
        for groups in [{.001: [6], .002: [3, 4, 5]},
                       {.001: [4, 5, 6], .002: [3], .003: [1, 2, 3]}]:
            result = analyze_bud_counts(observations(groups), list(groups))
            self.assertTrue(result.p_raw.isna().all())
            self.assertTrue(result.p_adjusted.isna().all())
            self.assertTrue((result.annotation == "n/a").all())

    def test_two_groups_match_pooled_student_t(self):
        groups = {.001: [6, 5, 6], .002: [1, 1, 3]}
        result = analyze_bud_counts(observations(groups), list(groups))
        reference = stats.ttest_ind(groups[.002], groups[.001], equal_var=True)
        self.assertAlmostEqual(result.t_comparison.iloc[0], reference.statistic)
        self.assertAlmostEqual(result.p_adjusted.iloc[0], reference.pvalue)

    def test_reference_is_thinnest_even_when_inputs_are_shuffled(self):
        groups = {.003: [4, 5], .001: [1, 2], .002: [2, 3]}
        data = observations(groups)
        expected = analyze_bud_counts(data, list(groups))
        observed = analyze_bud_counts(data.sample(frac=1, random_state=2), list(reversed(groups)))
        pd.testing.assert_frame_equal(expected, observed)
        self.assertTrue((observed.baseline_z == .001).all())

    def test_invalid_counts_and_duplicate_simulations_are_rejected(self):
        valid = observations({.001: [1, 2], .002: [3, 4]})
        for value in [-1, 0, 1.5, np.nan]:
            data = valid.astype({"N_buds": float}); data.loc[0, "N_buds"] = value
            with self.assertRaises(ValueError):
                analyze_bud_counts(data, [.001, .002])
        with self.assertRaises(ValueError):
            analyze_bud_counts(pd.concat([valid, valid.iloc[:1]]), [.001, .002])


if __name__ == "__main__":
    unittest.main()
