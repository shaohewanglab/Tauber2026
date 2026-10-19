#!/usr/bin/env python3
"""Decode the bundled blinded counts and match glands across the two intervals."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from bud_count_io import load_pairs
from bud_count_stats import prepare_intervals

FIGURE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = FIGURE_ROOT / 'source_data/washout-inputs.json'


def interval_fingerprint(intervals):
    """Hash matching, times, and integer endpoints independently of CSV formatting."""
    rows = []
    for row in intervals.sort_values(['subject_id', 'start_h']).itertuples():
        rows.append({'subject_id': row.subject_id, 'group': row.group, 'gland': str(row.gland),
                     'start_h': float(row.start_h), 'end_h': float(row.end_h),
                     'bud_start': int(row.bud_start), 'bud_end': int(row.bud_end)})
    return hashlib.sha256(json.dumps(rows, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def prepare_dataset(entry, figure_root=FIGURE_ROOT):
    sources = []
    for field in ('record', 'counts'):
        source = (figure_root / entry[field]).resolve()
        if not source.is_file():
            raise ValueError(f'Missing source file: {source}')
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        if digest != entry[field + '_sha256']:
            raise ValueError(f'Source checksum differs from the manifest: {entry[field]}')
        sources.append(source)
    data = load_pairs([sources], experiment=entry['experiment'])
    if len(data) != entry['expected_images'] or data.subject_id.nunique() != entry['expected_glands']:
        raise ValueError(f'Unexpected image or gland count for {entry["experiment"]}.')
    intervals = prepare_intervals(data, times=entry['times_h'])
    if interval_fingerprint(intervals) != entry['expected_interval_fingerprint']:
        raise ValueError(f'Matched gland counts differ from the validated figure data: {entry["experiment"]}.')
    # Portable provenance, relative to BB-94-washout rather than the local machine.
    data['record_file'] = entry['record']
    data['count_file'] = entry['counts']
    return data, intervals


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=DEFAULT_MANIFEST,
                        help='Input manifest; its source paths are relative to BB-94-washout.')
    parser.add_argument('--output-root', type=Path, default=FIGURE_ROOT / 'results',
                        help='Destination containing one subfolder per dataset (default: ../results).')
    args = parser.parse_args(argv)
    output_root = args.output_root.expanduser().resolve()
    for protected in (FIGURE_ROOT / 'source_data', FIGURE_ROOT / 'plotting_scripts'):
        if output_root == protected or protected in output_root.parents:
            parser.error('Output must not be inside source_data or plotting_scripts.')
    try:
        manifest = json.loads(args.manifest.read_text())
        # Validate every source before replacing generated results.
        prepared = [(entry, *prepare_dataset(entry)) for entry in manifest['datasets']]
        for entry, decoded, intervals in prepared:
            output = output_root / entry['experiment']
            output.mkdir(parents=True, exist_ok=True)
            decoded.to_csv(output / 'washout-bud-count-decoded.csv', index=False)
            intervals.to_csv(output / 'washout-bud-count-intervals.csv', index=False)
            record = {'experiment': entry['experiment'], 'times_h': entry['times_h'],
                      'image_counts': len(decoded), 'glands': decoded.subject_id.nunique(),
                      'source_paths_relative_to': 'BB-94-washout',
                      'inputs': [{'path': entry[field], 'sha256': entry[field + '_sha256']}
                                 for field in ('record', 'counts')],
                      'interval_fingerprint': interval_fingerprint(intervals),
                      'validation': 'Matches original validated gland identities, observation times, and integer endpoints.',
                      'timing_note': ('Recorded hours.' if set(decoded.time_unit) == {'hour'} else
                                      'Day labels converted to nominal hours (day x 24); exact elapsed hours are unavailable.')}
            (output / 'washout-bud-count-run.json').write_text(json.dumps(record, indent=2) + '\n')
            print(f'{entry["experiment"]}: {len(decoded)} images, {decoded.subject_id.nunique()} glands, {len(intervals)} matched interval ratios.')
    except (KeyError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
