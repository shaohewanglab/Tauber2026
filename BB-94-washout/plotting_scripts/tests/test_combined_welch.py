"""Independent checks of Welch variance, degrees of freedom and Holm family."""
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import unittest

import numpy as np
import pandas as pd
from scipy import stats

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


class WelchChecks(unittest.TestCase):
    def test_t_and_welch_df_from_group_moments(self):
        data = synthetic_data()
        result = comparisons(data)
        for row, stage in zip(result.itertuples(), ['first', 'second']):
            a = data.loc[(data.stage == stage) & (data.group == 'ctrl'), 'log2_ratio'].to_numpy()
            b = data.loc[(data.stage == stage) & (data.group == 'BB94'), 'log2_ratio'].to_numpy()
            va, vb = a.var(ddof=1)/len(a), b.var(ddof=1)/len(b)
            t = (b.mean()-a.mean())/math.sqrt(va+vb)
            df = (va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1))
            self.assertAlmostEqual(row.t_BB94_minus_control, t)
            self.assertAlmostEqual(row.df, df)
            self.assertAlmostEqual(row.p_raw, 2*stats.t.sf(abs(t),df))
            half = stats.t.ppf(.975,df)*math.sqrt(va+vb)
            self.assertAlmostEqual(row.ci95_unadjusted_low, b.mean()-a.mean()-half)
            self.assertAlmostEqual(row.ci95_unadjusted_high, b.mean()-a.mean()+half)

    def test_holm_family_contains_only_two_between_treatment_tests(self):
        result = comparisons(synthetic_data())
        self.assertEqual(len(result), 2)
        order = np.argsort(result.p_raw)
        expected = np.empty(2)
        expected[order[0]] = min(1, 2*result.p_raw.iloc[order[0]])
        expected[order[1]] = max(expected[order[0]], result.p_raw.iloc[order[1]])
        np.testing.assert_allclose(result.p_Holm, expected)
        shuffled = comparisons(synthetic_data().sample(frac=1,random_state=17))
        pd.testing.assert_frame_equal(result, shuffled)

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
