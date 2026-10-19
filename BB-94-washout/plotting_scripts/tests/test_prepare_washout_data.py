"""Guard explicit gland matching and agreement with the validated figure inputs."""
import json
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bud_count_stats import prepare_intervals
from prepare_washout_data import DEFAULT_MANIFEST, interval_fingerprint, prepare_dataset


class PreparationChecks(unittest.TestCase):
    def test_packaged_sources_reproduce_validated_intervals(self):
        manifest = json.loads(DEFAULT_MANIFEST.read_text())
        sizes = []
        for entry in manifest['datasets']:
            decoded, intervals = prepare_dataset(entry)
            sizes.append((len(decoded), len(intervals), intervals.subject_id.nunique()))
            self.assertEqual(interval_fingerprint(intervals), entry['expected_interval_fingerprint'])
            np.testing.assert_allclose(intervals.ratio, intervals.bud_end / intervals.bud_start)
            np.testing.assert_allclose(intervals.log2_ratio, np.log2(intervals.ratio))
        self.assertEqual(sizes, [(48, 32, 16), (78, 52, 26)])

    def test_shuffled_rows_match_by_identity_not_position(self):
        records = []
        for group, sequence in [('ctrl', [2, 10, 40]), ('BB94', [3, 6, 30])]:
            for time, count in zip([0, 24, 48], sequence):
                records.append({'experiment': 'test', 'group': group, 'gland': '1',
                                'time_h': time, 'bud_number': count})
        data = pd.DataFrame(records).sample(frac=1, random_state=73)
        intervals = prepare_intervals(data, times=(0, 24, 48))
        self.assertEqual(intervals.ratio.tolist(), [5, 2, 4, 5])
        self.assertEqual(intervals.subject_id.nunique(), 2)
        with self.assertRaisesRegex(ValueError, 'missing observation'):
            prepare_intervals(data.iloc[:-1], times=(0, 24, 48))

    def test_no_pseudocount_and_no_duplicate_gland_time(self):
        entry = json.loads(DEFAULT_MANIFEST.read_text())['datasets'][0]
        decoded, _ = prepare_dataset(entry)
        zero = decoded.copy()
        zero.loc[0, 'bud_number'] = 0
        with self.assertRaisesRegex(ValueError, 'positive'):
            prepare_intervals(zero, times=entry['times_h'])
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            prepare_intervals(pd.concat([decoded, decoded.iloc[[0]]]), times=entry['times_h'])


if __name__ == '__main__':
    unittest.main()
