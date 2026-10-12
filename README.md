# Tauber 2026 modeling: plot reproduction and model training

This repository contains the scripts and measurement tables for the Figure 4 and
Figure 5 experimental plots, plus training code and small test datasets for
3D-RCAN denoising and Cellpose-SAM segmentation of SMG live-imaging data.

## Quick start: Figure 4 and Figure 5 plots

The following environment and commands reproduce the experimental plots. For the
models, follow the separate [training instructions](#model-training-and-test-data).

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

3D-RCAN-training/                      Denoising configuration, modified inference script, and test data.
cellpose-SAM-training/                 Segmentation training notebooks, test script, and test data.

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

## Model training and test data

The training folders contain code, installation instructions, and small test
crops for the two models used to process the SMG live-imaging data. The full
training datasets and trained weights are deposited on Figshare in the
collection for Tauber et al., 2026. The included test data are crops from those
datasets.

| Folder and instructions | Model and purpose | Figshare items |
| --- | --- | --- |
| [3D-RCAN training](3D-RCAN-training/README.md) | **RCAN-denoise-SMG-BM**: denoising of both live-imaging channels | `3D-RCAN-training`, `trained-models` |
| [Cellpose-SAM training](cellpose-SAM-training/README.md) | **SMG-BM-v0907**: production epithelium segmentation; **SMG-BM-15V**: benchmarking | `cellpose-SAM-training`, `trained-models` |

Each folder's README covers installation, a quick test using the included data,
and reproduction of training. Follow the model-specific environments in the
[3D-RCAN installation instructions](3D-RCAN-training/rcan_installation.md) and
[Cellpose-SAM installation instructions](cellpose-SAM-training/cellpose_installation.md).
The root `requirements.txt` and `reproduce_plots.py` are for plotting only.

The modified [`3D-RCAN-training/apply.py`](3D-RCAN-training/apply.py) retains its
upstream **CC BY-NC 4.0** license, including its noncommercial restriction.
[`voxel_size_utils.py`](3D-RCAN-training/voxel_size_utils.py), the Cellpose-SAM
notebooks, and its test script are MIT-licensed. The bundled test data are covered
by CC BY 4.0; see [Licensing](LICENSE.md) for the scope and attribution.

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

## Plot source data and reproducibility

Nineteen original source tables are included byte-for-byte: eight Fig. 4
bud-count/key files, seven collagen IV line-scan CSVs, and four Fig. 5
bud-count/key files. Figure 5 ratios are rebuilt from the raw count/key pairs;
the scripts do not rely on copied intermediate ratio exports.

The plotting workflow uses these measurement tables and does not require image
data, trained models, modeling simulations, or manuscript files. Small imaging
test datasets are included separately in the two training folders.

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
| Original code, notebooks, software configuration, and documentation; exception below | [MIT](LICENSE-CODE) |
| Modified `3D-RCAN-training/apply.py` | [CC BY-NC 4.0](3D-RCAN-training/LICENSE-APPLY), noncommercial use only |
| Source and test data, derived data, plots, and statistical reports | [CC BY 4.0](LICENSE-DATA) |

See [License scope and attribution](LICENSE.md) for the files covered by each
license and how to credit the data.
