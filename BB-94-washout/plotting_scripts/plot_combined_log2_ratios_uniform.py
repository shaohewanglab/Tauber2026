#!/usr/bin/env python3
"""Combined plot with two Holm-corrected Welch t comparisons."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import plot_style as style

from plot_combined_log2_ratios import HERE, RESULTS, CATEGORIES, load_data
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

STEM = 'Fig5C-washout-log2-ratios-uniform-stars'


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
        result = stats.ttest_ind(treated, control, alternative='two-sided', equal_var=False)
        if not np.isfinite(result.statistic) or not np.isfinite(result.pvalue):
            raise ValueError('Welch inference requires a positive standard error.')
        ci = result.confidence_interval(confidence_level=.95)
        rows.append({'interval':label, 'test':'Two-sided Welch t-test', 'n_control':len(control), 'n_BB94':len(treated),
                     'mean_control':control.mean(), 'mean_BB94':treated.mean(),
                     'sd_control':control.std(ddof=1), 'sd_BB94':treated.std(ddof=1),
                     'mean_difference_BB94_minus_control':treated.mean()-control.mean(),
                     't_BB94_minus_control':result.statistic, 'df':result.df,
                     'ci95_unadjusted_low':ci.low, 'ci95_unadjusted_high':ci.high,
                     'p_raw':result.pvalue})
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
    rc = style.rc(args.font)
    size_key = 'fig5-washout-symbols' if dataset_markers else 'fig5-washout-uniform'
    import seaborn as sns
    with matplotlib.rc_context(rc):
        fig = plt.figure(figsize=style.dimensions(size_key), dpi=300)
        ax = fig.add_axes([0.12,0.1,0.85,0.78])
        sns.swarmplot(data=data, x='category', y='log2_ratio', order=CATEGORIES,
                      color='black', marker='o',
                      size=1.3 if dataset_markers else style.COMPACT_BUD_RATIO_CIRCLE_DIAMETER,
                      linewidth=0 if dataset_markers else style.COMPACT_CIRCLE_EDGE_WIDTH,
                      ax=ax, warn_thresh=0)
        ax.set_ylim(-0.5,4)
        ax.set_yticks([0,2,4])
        ax.set_ylabel(r'$\mathregular{Log_2}$(bud ratio)')
        ax.set_xlabel('')
        ax.set_xticks(range(4), labels=['Control d1/d0', 'BB-94 d1/d0',
                                      'Control d2/d1', 'BB-94 d2/d1'], rotation=45,
                      ha='right', rotation_mode='anchor')
        # Add breathing room between the outer category swarms and the frame.
        ax.set_xlim(-0.65, 3.65)
        fig.canvas.draw()
        collections = list(ax.collections)
        assert len(collections) == 4
        for index, (category, collection) in enumerate(zip(CATEGORIES, collections)):
            values = data.loc[data.category == category, 'log2_ratio']
            np.testing.assert_allclose(collection.get_offsets()[:,1], values)
            color = style.palettes()['washout']['control' if category.startswith('ctrl_') else 'BB94']
            if dataset_markers:
                collection.set_facecolor(color)
                collection.set_edgecolor(color)
                paths = []
                for marker in data.loc[data.category == category, 'marker']:
                    marker_style = MarkerStyle(marker)
                    paths.append(marker_style.get_path().transformed(marker_style.get_transform()))
                collection.set_paths(paths)
            else:
                style.compact_circle_style(collection, color)
            summary_color = 'black' if category.startswith('ctrl_') else 'gray'
            ax.errorbar(index, values.mean(), yerr=values.std(ddof=1), fmt='_', color=summary_color,
                        markersize=3, elinewidth=style.ERROR_WIDTH, capsize=1,
                        capthick=style.ERROR_WIDTH, markeredgewidth=style.MEAN_WIDTH,
                        zorder=4 if dataset_markers else style.COMPACT_SUMMARY_ZORDER)
        bracket_y = max(3.3, float(data.log2_ratio.max()) + 0.25)
        assert bracket_y + 0.3 < 4
        for pair, annotation in zip([(0,1), (2,3)], results.annotation):
            left, right = pair
            ax.plot([left,left,right,right], [bracket_y-0.09,bracket_y,bracket_y,bracket_y-0.09], color='black', lw=style.FRAME_WIDTH)
            ax.text((left+right)/2, bracket_y+0.02, annotation, ha='center', va='bottom', fontsize=style.ANNOTATION_SIZE)
        style.compact_layout(fig, ax, size_key, vector_mathtext=True)
        frame_inches = (ax.get_position().size * fig.get_size_inches()).tolist()
        style.save_native_figure(fig, output / stem, dpi=450)
        plt.close(fig)
    results.to_csv(output/(stem+'-statistics.csv'), index=False)
    data.to_csv(output / 'Fig5C-washout-log2-ratios-data.csv', index=False)
    data.groupby(['dataset', 'category']).log2_ratio.agg(n='size', mean='mean', sd='std').reset_index().to_csv(
        output / 'Fig5C-washout-log2-ratios-by-dataset.csv', index=False)
    data.groupby('category').log2_ratio.agg(n='size', mean='mean', sd='std').reindex(CATEGORIES).to_csv(
        output / 'Fig5C-washout-log2-ratios-summary.csv')
    record = {'inputs':provenance, 'glands_total':42, 'glands_per_group':21,
              'method':'Two-sided Welch t-tests on per-gland log2 ratios; Holm correction across the two requested control-versus-BB94 comparisons.',
              'p_value_calculation':'SciPy ttest_ind(equal_var=False, alternative=two-sided), with Welch-Satterthwaite degrees of freedom. No equal-variance assumption.',
              'effect_size':'Difference in mean log2 ratios, BB-94 minus control. Positive values indicate larger BB-94 mean ratios on the log2 scale; reported 95% confidence intervals are pointwise, not multiplicity-adjusted.',
              'interpretation':'The test compares mean per-gland log2 ratios. Assumes independent groups and approximately normal errors or an adequate mean approximation; permits unequal variances.',
              'family_size':2,'alpha':0.05,'error_bars':'SD',
              'markers':('Filled circles: 210325-27; filled triangles: 210401-03; gray control, red BB-94.' if dataset_markers else 'Hollow circles with white interiors; gray control outlines, red BB-94 outlines; mean/SD drawn above points.'),
              'marker_diameter_points':1.3 if dataset_markers else style.COMPACT_BUD_RATIO_CIRCLE_DIAMETER,
              'marker_edge_width_points':0 if dataset_markers else style.COMPACT_CIRCLE_EDGE_WIDTH,
              'dataset_legend_displayed':False, 'font':args.font,
              'panel':'Figure 5C', 'colors':style.palettes()['washout'],
              'summary_colors':{'control':'black', 'BB94':'gray'},
              'x_tick_rotation_degrees':45, 'plotting_frame_inches':frame_inches,
              'x_axis_limits':list(ax.get_xlim()),
              'independence':'Different treatment groups are treated as independent. Each test uses one observation per gland; repeated use across intervals is covered by the Holm comparison family, not treated as extra glands within either test.',
              'pooling':'Glands from both datasets are pooled with equal gland weights. Dataset and plate are not fitted as covariates.',
              'alignment':'22h/0h with day1/day0; 46h/22h with day2/day1. Original ratios unchanged; no rate normalization.',
              'annotation_thresholds':{'*':'p_Holm < 0.05','**':'p_Holm < 0.01','***':'p_Holm < 0.001','****':'p_Holm < 0.0001','ns':'p_Holm >= 0.05'},
              'results':results.to_dict(orient='records')}
    (output/(stem+'-record.json')).write_text(json.dumps(record,indent=2)+'\n')
    paragraphs = ['COMBINED WASHOUT: CONTROL VERSUS BB-94',
                  '21 glands per treatment, comprising 8 from 210325-27 and 13 from 210401-03. Each point is one gland. Error bar is SD.',
                  'Two-sided Welch t-tests on per-gland log2 ratios, with Holm correction across these two requested comparisons. This is a two-comparison family; within-treatment interval comparisons are not included in this figure or analysis.',
                  'Welch t-tests compare mean log2 ratios using separate group variances and Welch-Satterthwaite degrees of freedom. Independent groups and approximately normal errors or an adequate mean approximation are assumed. Confidence intervals for mean differences are pointwise 95%, not multiplicity-adjusted.',
                  'Treatment groups are independent; gland matching across time is preserved in calculation of each ratio. No gland contributes more than once to either individual test.',
                  'Data are pooled by experimental stage: 22h/0h corresponds to day1/day0; 46h/22h corresponds to day2/day1. Ratios are not rate-normalized. Day annotations do not establish precise elapsed hours.',
                  'This is a pooled gland-level comparison; dataset and plate are not covariates in the Welch tests.']
    for row in results.itertuples():
        paragraphs.append(f'{row.interval}: {row.annotation}\n'
                          f'Control mean log2 ratio = {row.mean_control:.6f}, SD = {row.sd_control:.6f}; BB-94 mean = {row.mean_BB94:.6f}, SD = {row.sd_BB94:.6f}.\n'
                          f'Mean difference (BB-94 minus control) = {row.mean_difference_BB94_minus_control:.6f}; pointwise 95% CI [{row.ci95_unadjusted_low:.6f}, {row.ci95_unadjusted_high:.6f}].\n'
                          f'Welch t = {row.t_BB94_minus_control:.6f}; df = {row.df:.6f}.\n'
                          f'Raw p = {row.p_raw:.10g}; Holm-adjusted p = {row.p_Holm:.10g}.')
    paragraphs.append('Asterisks use Holm-adjusted p: * < .05; ** < .01; *** < .001; **** < .0001.')
    (output/(stem+'-statistics.txt')).write_text('\n\n'.join(paragraphs)+'\n')
    print(results[['interval','p_raw','p_Holm','annotation']].to_string(index=False))
    print('Saved '+stem+' as PDF, SVG, PNG, plus statistics CSV/TXT and method record.')


if __name__ == '__main__':
    main()
