# Tauber 2026 modeling: Figure 4 and Figure 5 plot reproduction

This folder contains the scripts and all measurement tables needed to reproduce
the current experimental plots.

## Quick start

Use Python 3.11 or newer (tested with Python 3.12.11). From this folder:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce_plots.py
```

The existing scientific Python environment can be used without creating a venv
if it already has the listed packages. The requirements file pins tested versions.
On Windows, activate `.venv` with `.venv\Scripts\activate` instead.

Arial is required to reproduce the manuscript typography. No font is bundled.
For a system without Arial:

```bash
python reproduce_plots.py --font "DejaVu Sans"
```

That changes typography only. All numerical checks remain the same. All paths
are resolved within this folder, so it can be moved or cloned elsewhere.
The workflow replaces generated results; it does not modify the source data.
It does not require Jupyter, the manuscript directory, or a network connection
after Python dependencies are installed.

## Contents

```text
Fig4-BB94-titration/
  source_data/                         Original blinded count/key tables and collagen IV CSVs.
  plotting_scripts/                    Two plotting CLIs, input manifests, shared modules/tests.
  results/bud-count/                   Bud-count plot, individual gland data and statistics.
  results/collagen-IV/                 Collagen IV plot, gland data and statistics.

Fig5-washout/
  source_data/                         Original count/key pairs for 210325-27 and 210401-03.
  plotting_scripts/                    Decode/ratio preparation and both combined plot variants.
  results/                             Rebuilt per-dataset ratios and combined plots/reports.

source-data-manifest.json              SHA-256 hashes and original relative source names.
validation/expected-results.json       Numerical reference from accepted analyses.
validation/last-verification.json      Latest reproduction verification result.
reproduce_plots.py                     Rebuild all four current plots from source tables.
verify_reproduction.py                 Verify input hashes and accepted numerical results.
```

Each plot is exported as PDF, editable-text SVG and PNG. Readable statistical
reports and full-precision CSV results accompany the plots. Existing output
names are preserved. Fig. 5 has a uniform-symbol version and a version with
dataset-specific symbols, without a dataset-name legend.

## Statistical methods

Both figures use two-sided Mann-Whitney U tests with Holm correction. Each point
is one gland. Horizontal bars indicate means; error bars are SD. Mann-Whitney
compares ranks/distributions rather than directly testing means. A pure median
interpretation additionally requires similar distribution shapes.

Fig. 4: compare each of the six nonzero doses with DMSO (0). Holm correction is
applied separately across those six comparisons for each plotted outcome.
No treatment-versus-treatment comparisons enter the current correction family.

Bud counts use `log2(day1/day0)`, pairing each gland across days first. There are
110 glands from three cohorts (33, 49 and 28). The 0.005 and 0.05 uM groups have
15 glands each; other doses have 16. Two filename-marked placeholder glands are
excluded explicitly. Future-washout day2 observations are outside this endpoint.

Collagen IV uses the mean of each gland's retained line-scan AUC measurements:
716 scans, 49 glands, seven glands per dose, 12-15 scans per gland. Divide the
gland means by the mean of the seven DMSO gland means, making control mean 1.
These source tables are measurements from raw images; no denoising is applied
by the plotting workflow. The statistical unit is the gland, not the line scan.

Fig. 5: compare control with BB-94 for each of two successive intervals; apply
Holm correction across those two comparisons. There are 21 glands per treatment
(8 + 13 from two datasets), with matching across time retained in each ratio.
Recorded 0/22/46 hours and annotated day0/day1/day2 are aligned by experimental
stage; ratios are not converted to rates. Symbols represent datasets when used.

Exact Mann-Whitney p-values are used for small untied groups (collagen IV);
bud-count comparisons use asymptotic calculations with tie and continuity
corrections. Stars use adjusted p: `*` < .05; `**` < .01; `***` < .001; `****` < .0001.
All analyses treat distinct glands as independent across treatment groups.

## Source data and reproducibility

Nineteen original source tables are included byte-for-byte: eight Fig. 4
bud-count/key files, seven collagen IV line-scan CSVs, and four Fig. 5
bud-count/key files. Figure 5 ratios are rebuilt from the raw count/key pairs;
the scripts do not rely on copied intermediate ratio exports.

These are the source measurements for plotting. Raw microscope images, image
segmentation/training data, modeling simulations and manuscript files are not
needed for these plots and are outside this package.

The verifier checks all 19 source hashes and compares regenerated gland values,
summaries and statistics to the accepted analysis. It runs automatically after
reproduction, or separately with:

```bash
python verify_reproduction.py
```

For unit tests:

```bash
python -m unittest discover -s Fig4-BB94-titration/plotting_scripts -p "test_*.py"
python -m unittest discover -s Fig5-washout/plotting_scripts/tests -p "test_*.py"
```

Paths in generated provenance records are stored relative to this package when
running `reproduce_plots.py`. Individual CLIs may report absolute paths for their
current execution location; these are audit information, not external inputs.

See the [Figure 4 CLI instructions](Fig4-BB94-titration/plotting_scripts/README-titration-CLI.txt)
and [Figure 5 instructions](Fig5-washout/README.txt) for individual commands
and file formats.

## Method references

- [SciPy: Mann–Whitney U test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html)
- [statsmodels: multiple-testing corrections](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html)

## License

| Material | License |
| --- | --- |
| Code, software configuration, and associated documentation | [MIT](LICENSE-CODE) |
| Source data, derived data, plots, and statistical reports | [CC BY 4.0](LICENSE-DATA) |

See [License scope and attribution](LICENSE.md) for the files covered by each
license and how to credit the data.
