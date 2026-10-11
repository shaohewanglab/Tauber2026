#!/usr/bin/env python3
"""Combined plot with two Holm-corrected Mann-Whitney U comparisons."""
import argparse
import json
from pathlib import Path

from plot_combined_log2_ratios import HERE, RESULTS, CATEGORIES, LABELS, load_data
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

STEM = 'combined-washout-log2-ratios-uniform-stars'


def stars(p):
    if not np.isfinite(p):
        raise ValueError('Cannot annotate a non-finite p-value.')
    for threshold, label in [(0.0001, '****'), (0.001, '***'), (0.01, '**'), (0.05, '*')]:
        if p < threshold:
            return label
    return 'ns'


def comparisons(data):
    required = {'subject_id', 'stage', 'group', 'log2_ratio'}
    if not required.issubset(data.columns):
        raise ValueError('Missing comparison columns: ' + ', '.join(sorted(required - set(data.columns))))
    if data.duplicated(['subject_id', 'stage']).any():
        raise ValueError('Each gland must contribute only once per interval.')
    if set(data.group) != {'ctrl', 'BB94'} or set(data.stage) != {'first', 'second'}:
        raise ValueError('Expected control and BB94 for both intervals.')
    if (data.groupby('subject_id').group.nunique() != 1).any():
        raise ValueError('A gland cannot belong to both treatment groups.')
    rows = []
    for stage, label in [('first', 'd1/d0'), ('second', 'd2/d1')]:
        control = data.loc[(data.stage == stage) & (data.group == 'ctrl'), 'log2_ratio'].to_numpy()
        treated = data.loc[(data.stage == stage) & (data.group == 'BB94'), 'log2_ratio'].to_numpy()
        if min(len(control), len(treated)) < 2 or not np.isfinite(np.r_[control, treated]).all():
            raise ValueError('Each group needs at least two finite gland-level ratios.')
        tied = len(np.unique(np.r_[control, treated])) < len(control) + len(treated)
        method_used = 'exact' if min(len(control), len(treated)) <= 8 and not tied else 'asymptotic'
        result = stats.mannwhitneyu(treated, control, alternative='two-sided', method='auto', use_continuity=True)
        rows.append({'interval':label, 'test':'Two-sided Mann-Whitney U test', 'n_control':len(control), 'n_BB94':len(treated),
                     'mean_control':control.mean(), 'mean_BB94':treated.mean(),
                     'sd_control':control.std(ddof=1), 'sd_BB94':treated.std(ddof=1),
                     'descriptive_mean_difference_BB94_minus_control':treated.mean()-control.mean(),
                     'U_BB94':result.statistic,
                     'rank_biserial_BB94_minus_control':2 * result.statistic / (len(control) * len(treated)) - 1,
                     'method_used':method_used, 'ties_present':tied,
                     'continuity_correction':method_used == 'asymptotic', 'p_raw':result.pvalue})
    results = pd.DataFrame(rows)
    results['p_Holm'] = multipletests(results.p_raw, method='holm', alpha=0.05)[1]
    results['annotation'] = results.p_Holm.map(stars)
    return results


def main(dataset_markers=False, stem=STEM, argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=RESULTS / 'combined',
                        help='Destination for the plot and report files (default: ../results/combined).')
    parser.add_argument('--font', default='Arial', help='Installed font family; use DejaVu Sans for a portable alternative.')
    args = parser.parse_args(argv)
    output = args.output_dir.expanduser().resolve()
    # Output may be relocated, but should never replace bundled inputs or scripts.
    for protected in (HERE, HERE.parent / 'source_data'):
        if output == protected or protected in output.parents:
            parser.error('Output must not be inside plotting_scripts or source_data.')
    output.mkdir(parents=True, exist_ok=True)
    data, provenance = load_data()
    results = comparisons(data)
    rc = {'font.family':args.font, 'font.size':7, 'svg.fonttype':'none', 'pdf.fonttype':42, 'axes.linewidth':0.8}
    import seaborn as sns
    with matplotlib.rc_context(rc):
        fig = plt.figure(figsize=(1.9,1.65), dpi=300)
        ax = fig.add_axes([0.12,0.1,0.85,0.78])
        sns.swarmplot(data=data, x='category', y='log2_ratio', order=CATEGORIES,
                      color='black', marker='o', size=2.7, ax=ax, warn_thresh=0)
        ax.set_ylim(-0.5,4)
        ax.set_yticks([0,2,4])
        ax.set_ylabel(r'$\mathregular{Log_2}$(bud ratio)')
        ax.set_xlabel('')
        ax.set_xticks(range(4), labels=LABELS, rotation=45, ha='right')
        fig.canvas.draw()
        collections = list(ax.collections)
        assert len(collections) == 4
        for index, (category, collection) in enumerate(zip(CATEGORIES, collections)):
            values = data.loc[data.category == category, 'log2_ratio']
            np.testing.assert_allclose(collection.get_offsets()[:,1], values)
            color = 'black' if category.startswith('ctrl_') else 'red'
            collection.set_facecolor(color)
            collection.set_edgecolor(color)
            if dataset_markers:
                paths = []
                for marker in data.loc[data.category == category, 'marker']:
                    style = MarkerStyle(marker)
                    paths.append(style.get_path().transformed(style.get_transform()))
                collection.set_paths(paths)
            ax.errorbar(index, values.mean(), yerr=values.std(ddof=1), fmt='_', color='gray',
                        markersize=5, elinewidth=0.7, capsize=2, capthick=0.7, zorder=4)
        bracket_y = max(3.3, float(data.log2_ratio.max()) + 0.25)
        assert bracket_y + 0.3 < 4
        for pair, annotation in zip([(0,1), (2,3)], results.annotation):
            left, right = pair
            ax.plot([left,left,right,right], [bracket_y-0.09,bracket_y,bracket_y,bracket_y-0.09], color='black', lw=0.7)
            ax.text((left+right)/2, bracket_y+0.02, annotation, ha='center', va='bottom', fontsize=9)
        for suffix in ['svg','pdf','png']:
            fig.savefig(output / (stem+'.'+suffix), dpi=450, bbox_inches='tight', pad_inches=0.04)
        plt.close(fig)
    results.to_csv(output/(stem+'-statistics.csv'), index=False)
    data.to_csv(output / 'combined-washout-log2-ratios-data.csv', index=False)
    data.groupby(['dataset', 'category']).log2_ratio.agg(n='size', mean='mean', sd='std').reset_index().to_csv(
        output / 'combined-washout-log2-ratios-by-dataset.csv', index=False)
    data.groupby('category').log2_ratio.agg(n='size', mean='mean', sd='std').reindex(CATEGORIES).to_csv(
        output / 'combined-washout-log2-ratios-summary.csv')
    record = {'inputs':provenance, 'glands_total':42, 'glands_per_group':21,
              'method':'Two-sided Mann-Whitney U tests on per-gland log2 ratios; Holm correction across the two requested control-versus-BB94 comparisons.',
              'p_value_calculation':'SciPy method=auto, use_continuity=True. Both current comparisons use the asymptotic method with tie and continuity corrections.',
              'effect_size':'Rank-biserial correlation = 2*U_BB94/(n_BB94*n_control)-1. Positive values indicate larger BB-94 ratios; ties contribute half to U.',
              'interpretation':'The test compares the two distributions. Mean differences are descriptive; this is not specifically a test of means or medians.',
              'family_size':2,'alpha':0.05,'error_bars':'SD',
              'markers':('Circles: 210325-27; triangles: 210401-03; black control, red BB-94.' if dataset_markers else 'Uniform circles; black control, red BB-94.'),
              'dataset_legend_displayed':False, 'font':args.font,
              'independence':'Different treatment groups are treated as independent. Each test uses one observation per gland; repeated use across intervals is covered by the Holm comparison family, not treated as extra glands within either test.',
              'pooling':'Glands from both datasets are pooled with equal gland weights. Dataset and plate are not fitted as covariates.',
              'alignment':'22h/0h with day1/day0; 46h/22h with day2/day1. Original ratios unchanged; no rate normalization.',
              'annotation_thresholds':{'*':'p_Holm < 0.05','**':'p_Holm < 0.01','***':'p_Holm < 0.001','****':'p_Holm < 0.0001','ns':'p_Holm >= 0.05'},
              'results':results.to_dict(orient='records')}
    (output/(stem+'-record.json')).write_text(json.dumps(record,indent=2)+'\n')
    paragraphs = ['COMBINED WASHOUT: CONTROL VERSUS BB-94',
                  '21 glands per treatment, comprising 8 from 210325-27 and 13 from 210401-03. Each point is one gland. Error bar is SD.',
                  'Two-sided Mann-Whitney U tests on per-gland log2 ratios, with Holm correction across these two requested comparisons. This is a two-comparison family; within-treatment interval comparisons are not included in this figure or analysis.',
                  'SciPy method=auto chooses the asymptotic method for both current comparisons; tied ranks and continuity correction are included. The test compares distributions rather than specifically testing their means or medians. Means and SDs describe the plotted data.',
                  'Treatment groups are independent; gland matching across time is preserved in calculation of each ratio. No gland contributes more than once to either individual test.',
                  'Data are pooled by experimental stage: 22h/0h corresponds to day1/day0; 46h/22h corresponds to day2/day1. Ratios are not rate-normalized. Day annotations do not establish precise elapsed hours.',
                  'This is a pooled gland-level comparison; dataset and plate are not covariates in the Mann-Whitney tests.']
    for row in results.itertuples():
        paragraphs.append(f'{row.interval}: {row.annotation}\n'
                          f'Control mean log2 ratio = {row.mean_control:.6f}, SD = {row.sd_control:.6f}; BB-94 mean = {row.mean_BB94:.6f}, SD = {row.sd_BB94:.6f}.\n'
                          f'Descriptive mean difference (BB-94 minus control) = {row.descriptive_mean_difference_BB94_minus_control:.6f}.\n'
                          f'U (BB-94 as first sample) = {row.U_BB94:.6f}; rank-biserial correlation = {row.rank_biserial_BB94_minus_control:.6f}.\n'
                          f'Raw p = {row.p_raw:.10g}; Holm-adjusted p = {row.p_Holm:.10g}; method = {row.method_used}.')
    paragraphs.append('Rank-biserial correlation = 2*U_BB94/(n_BB94*n_control)-1; positive values indicate larger BB-94 ratios and negative values indicate smaller BB-94 ratios. Ties contribute half to U.')
    paragraphs.append('Asterisks use Holm-adjusted p: * < .05; ** < .01; *** < .001; **** < .0001.')
    (output/(stem+'-statistics.txt')).write_text('\n\n'.join(paragraphs)+'\n')
    print(results[['interval','p_raw','p_Holm','annotation']].to_string(index=False))
    print('Saved '+stem+' as PDF, SVG, PNG, plus statistics CSV/TXT and method record.')


if __name__ == '__main__':
    main()
