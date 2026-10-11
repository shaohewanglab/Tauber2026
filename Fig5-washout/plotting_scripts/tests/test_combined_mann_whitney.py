"""Independent checks of tied ranks, small-sample exact p, and Holm family."""
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import unittest

import numpy as np
import pandas as pd

from plot_combined_log2_ratios_uniform import comparisons, stars


def synthetic_data():
    rows = []
    values = {'first': {'ctrl': [1, 2, 3], 'BB94': [2, 2, 5]},
              'second': {'ctrl': [3, 4, 5], 'BB94': [0, 1, 2]}}
    for stage, groups in values.items():
        for group, measurements in groups.items():
            for i, value in enumerate(measurements):
                rows.append({'subject_id': f'{group}-{i}', 'stage': stage,
                             'group': group, 'log2_ratio': value})
    return pd.DataFrame(rows)


class MannWhitneyChecks(unittest.TestCase):
    def test_tied_u_asymptotic_p_and_effect_by_hand(self):
        result = comparisons(synthetic_data()).iloc[0]
        # U = 1.5 + 1.5 + 3 = 6; tie block has size 3.
        variance = (3 * 3 / 12) * (7 - (3**3 - 3) / (6 * 5))
        z = (6 - 4.5 - 0.5) / math.sqrt(variance)
        p = math.erfc(z / math.sqrt(2))
        self.assertEqual(result.U_BB94, 6)
        self.assertAlmostEqual(result.p_raw, p, places=14)
        self.assertAlmostEqual(result.rank_biserial_BB94_minus_control, 1 / 3)
        self.assertEqual(result.method_used, 'asymptotic')
        self.assertTrue(result.ties_present)
        self.assertTrue(result.continuity_correction)

    def test_exact_separation_and_holm_family_by_hand(self):
        result = comparisons(synthetic_data())
        # Complete separation: 2 out of C(6,3)=20 assignments are as extreme.
        self.assertEqual(result.iloc[1].U_BB94, 0)
        self.assertAlmostEqual(result.iloc[1].p_raw, 2 / math.comb(6, 3))
        self.assertEqual(result.iloc[1].rank_biserial_BB94_minus_control, -1)
        self.assertEqual(result.iloc[1].method_used, 'exact')
        self.assertFalse(result.iloc[1].continuity_correction)
        self.assertAlmostEqual(result.iloc[1].p_Holm, 0.2)
        self.assertAlmostEqual(result.iloc[0].p_Holm, result.iloc[0].p_raw)

    def test_repeated_gland_is_not_an_extra_observation(self):
        data = synthetic_data()
        self.assertEqual(comparisons(data).n_control.tolist(), [3, 3])
        with self.assertRaisesRegex(ValueError, 'only once'):
            comparisons(pd.concat([data, data.iloc[[0]]], ignore_index=True))

    def test_nonfinite_and_cross_group_subjects_are_rejected(self):
        data = synthetic_data()
        data.loc[0, 'log2_ratio'] = np.nan
        with self.assertRaisesRegex(ValueError, 'finite'):
            comparisons(data)
        data = synthetic_data()
        data.loc[(data.stage == 'second') & (data.subject_id == 'BB94-0'), 'subject_id'] = 'ctrl-0'
        with self.assertRaisesRegex(ValueError, 'only once'):
            comparisons(data)

    def test_thresholds_are_strict(self):
        self.assertEqual(stars(0.05), 'ns')
        self.assertEqual(stars(0.01), '*')
        self.assertEqual(stars(0.001), '**')
        self.assertEqual(stars(0.0001), '***')
        self.assertEqual(stars(0.00001), '****')


if __name__ == '__main__':
    unittest.main()
