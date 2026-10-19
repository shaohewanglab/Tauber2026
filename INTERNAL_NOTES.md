# Internal manuscript data notes

## Figure 4E model bud counts - complete, 2026-10-18

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
