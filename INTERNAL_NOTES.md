# Internal manuscript data notes

## Figure 4A/E final inputs, 2026-11-10

This supersedes the October Figure 4 model inputs and statistics below. The source is OneDrive `ModelDataV2`: the selected 18 `thickness_paper/*_exp_V001-V003_k12.npz` archives and `bud_counts_summary_cusps.csv`. The current repository CSV is copied unchanged as `model-morphometry/ModelData/bud_counts_summary.csv`; its SHA256 is `f35534c76223fe166d4bc4a7eb625718bd3fbedcdc18ed7b222607da37fbb2e0`. The October table is retained under `model-morphometry/validation/history`.

| Effective z | Independent counts | Mean | Sample SD |
| --- | --- | ---: | ---: |
| 0.001 | 6, 7, 6 | 6.3333333333 | 0.5773502692 |
| 0.002 | 5, 6, 6 | 5.6666666667 | 0.5773502692 |
| 0.003 | 5, 5, 6 | 5.3333333333 | 0.5773502692 |
| 0.005 | 3, 4, 4 | 3.6666666667 | 0.5773502692 |
| 0.008 | 1, 1, 1 | 1 | 0 |
| 0.01 | 1, 1, 1 | 1 | 0 |

Shaohe requested mean +/- sample SD of all three independent simulations in 4A. Baseline each run separately, then summarize the change in acircularity. The metric, open-edge perimeter, area closure, 0.01 saved-frame interval and 0-1.5 time range are unchanged. The six selected effective z values and their colors match across 4A, its showcase shapes and 4E. Experimental 4B colors remain blue/green/red; their palette associations are updated to the second/fourth/sixth selected thicknesses (0.002/0.005/0.01), without inferring a new parameter fit.

The current five Dunnett comparisons versus z=0.001 give ns, ns, ****, ****, ****. Pooled within-group variance is 2/9, residual df=12. Adjusted P values are 0.3356949191, 0.0846301915, 0.00006922114055, 0.00000004304914446, 0.00000004304914446. The existing comparison method and assumptions are unchanged. The independent-run Spearman trend is rho=-0.9451620167, raw permutation P=0.000001, Holm P=0.000003 across the separate E/F/G trend family. Captions need the new baseline z and 4A mean/sample-SD description.

## Historical Figure 4E model bud counts, 2026-10-18

Justin's updated `ModelData/bud_counts_summary.csv` now contains all **18 intended
independent simulations**, three at each of six thicknesses. Counts are measured
at **tau/tau_g = 1.5**, as confirmed by Shaohe. The shared CSV is copied unchanged
into `model-morphometry/ModelData/bud_counts_summary.csv`.

| z | Counts (Nbud1, Nbud2, Nbud3) | Mean | Sample SD |
| --- | --- | ---: | ---: |
| 0.0005 | 6, 6, 6 | 6 | 0 |
| 0.001 | 5, 6, 5 | 5.3333333333 | 0.5773502692 |
| 0.003 | 4, 3, 4 | 3.6666666667 | 0.5773502692 |
| 0.005 | 1, 1, 3 | 1.6666666667 | 1.1547005384 |
| 0.008 | 1, 1, 1 | 1 | 0 |
| 0.01 | 1, 1, 1 | 1 | 0 |

Input SHA256: `e5e200abfdaca047712725192b0e46d552221117da84964b8f947dc32d793e97`. Both source manifests record this digest. All supplied
means and sample SDs are checked against the three raw counts before plotting.
The CSV column identities Nbud1-3 are retained; no simulation identifiers or
observations are invented.

Five two-sided, single-step Dunnett comparisons use z = 0.0005 as baseline.
All comparisons are available. Adjusted P values in ascending comparison-z order
are 0.5097920842, 0.001418209002, 0.000003913651762, 0.000000841973981
and 0.000000841973981; labels are ns, **, ****, ****, ****. Pooled within-group
variance is 1/3 with residual df = 12. Three simulations per condition cannot
establish the approximate normal/common-variance assumption; reports state it.
An unbudded rounded model shape counts as 1. No data are imputed.

The current 4E canvas remains 1.32 x 0.92 inches, with outlined Arial/Arev text,
Tol Bright hollow circles, black mean/SD drawn above observations, and the same
compact axis layout. Figure 4 uses only Paul Tol Bright; alternative color
exports are retired. The incomplete-data placeholder in the manuscript should
be replaced when Shaohe reaches that edit in the ongoing manuscript review;
this plot update does not change the manuscript review or its browser choices.

## History

On 17 October 2026 there were 14/18 counts: n=1 at z=0.0005 and z=0.001, and
n=3 at the four thicker values. Two additional runs were pending at each of the
two thinnest conditions. Their SDs and all five tests were unavailable (`n/a`).
The updated shared table received on 18 October completes those four counts;
no missing values were imputed at either stage.

## Statistical-method review, 2026-10-18

Shaohe requested retaining each-thickness-versus-thinnest annotations on Fig4E,
while adding a Spearman permutation trend test to the caption for 4E, 4F and 4G.
He confirmed that an unbudded model shape counts as one, not zero. On 18 October
he selected Dunnett for 4E, 4F and 4G. This supersedes the earlier rank-test and
shifted-Poisson method review. Dunnett families are separate for each panel;
there is no extra Holm correction. Figure 4 uses only Paul Tol Bright. The
historical exact-rank results are retained under validation/history, rather than
presented as the current inference. Experimental gland endpoints are unchanged.

All three requested Spearman trends have been implemented with 999,999
permutations, preserving ties, a two-sided absolute-rho statistic, and seed
20261018. Correlations are -0.938567 (E), -0.710800 (F), +0.886048 (G); the
Holm-adjusted trend P values are all 0.000003. See the root validation report.

## Figure 5C method selection, 2026-10-18

Shaohe selected two-sided Welch t-tests on per-gland log2 ratios, separately
comparing control with BB-94 during treatment and after washout. Holm correction
covers these two dependent interval comparisons. The same glands are followed
across time but each interval test uses independent treatment groups, one ratio
per gland. Both annotations remain ****. The original gland data are unchanged.
