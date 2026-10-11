#!/usr/bin/env python3
"""Load validated gland-level data for the two current combined Fig. 5 plots."""
import hashlib
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'bud-count-matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

HERE = Path(__file__).resolve().parent
FIGURE_ROOT = HERE.parent
RESULTS = FIGURE_ROOT / 'results'
SOURCES = [
    ('210325-27', RESULTS / '210325-washout/washout-bud-count-intervals.csv', [(0, 22), (22, 46)], 'o'),
    ('210401-03', RESULTS / '210401-washout/washout-bud-count-intervals.csv', [(0, 24), (24, 48)], '^'),
]
CATEGORIES = ['ctrl_first', 'BB94_first', 'ctrl_second', 'BB94_second']
LABELS = ['Control d1/d0', 'BB-94 d1/d0', 'Control washout d2/d1', 'BB-94 washout d2/d1']


def load_data():
    frames, provenance = [], []
    for label, source, times, marker in SOURCES:
        frame = pd.read_csv(source, dtype={'gland':str, 'subject_id':str})
        assert set(frame.group) == {'ctrl', 'BB94'}, source
        assert not frame.duplicated(['subject_id', 'start_h', 'end_h']).any(), source
        np.testing.assert_allclose(frame.ratio, frame.bud_end / frame.bud_start, rtol=1e-12)
        np.testing.assert_allclose(frame.log2_ratio, np.log2(frame.bud_end / frame.bud_start), rtol=1e-12, atol=1e-12)
        assert set(zip(frame.start_h, frame.end_h)) == set(times), source
        for subject, rows in frame.groupby('subject_id'):
            assert len(rows) == 2 and rows.group.nunique() == 1, subject
            rows = rows.sort_values('start_h')
            assert rows.iloc[0].bud_end == rows.iloc[1].bud_start, subject
        frame['dataset'] = label
        frame['marker'] = marker
        frame['stage'] = [dict(zip(times, ['first', 'second']))[(a,b)] for a,b in zip(frame.start_h, frame.end_h)]
        frame['category'] = frame.group + '_' + frame.stage
        frames.append(frame)
        provenance.append({'dataset':label, 'source':str(source.relative_to(FIGURE_ROOT)), 'sha256':hashlib.sha256(source.read_bytes()).hexdigest(), 'marker':marker,
                           'source_intervals_h':times, 'timing_note':'Recorded hours' if label=='210325-27' else 'Nominal day x 24; exact elapsed hours unavailable'})
    combined = pd.concat(frames, ignore_index=True)
    assert not combined.duplicated(['subject_id', 'stage']).any()
    assert combined.subject_id.nunique() == 42 and len(combined) == 84
    assert combined.groupby('category').size().reindex(CATEGORIES).tolist() == [21]*4
    return combined, provenance

