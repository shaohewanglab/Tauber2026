# Tauber 2026 modeling: plot reproduction and model training

This repository contains the source data and scripts for the Figure 2, Figure 4
and Figure 5 experimental and simulation plots, plus 3D-RCAN and Cellpose-SAM
training code and small imaging test datasets.

## Quick start: manuscript plots

Use Python 3.11 or newer (tested with Python 3.12.11):

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce_plots.py
```

This rebuilds **10 PDF/SVG/PNG plot sets**, numerical summaries and statistics,
then verifies bundled source hashes and accepted numerical results. All inputs
are local to this repository. The workflow requires Arial, TeX Live tools
`xelatex`, `xdvipdfmx`, `dvisvgm`, Poppler's `pdftoppm`, and Ghostscript's `gs`. The TeX packages are
`fontspec`, `amsmath`, `arevmath`, `pgf` and `standalone`.

Shared typography, line weights and layout are in `plot_style.py`; full canvas
and frame sizes are in `figure_sizes.json`; colors are in `plot_palettes.json`.
Text uses Arial and model/mesh mathematical symbols use Justin's Arev font.
Both PDF and SVG exports convert all text to vector outlines. No fonts need to
be installed to view or place these exports in Illustrator. Font installation
is needed only when regenerating plots; edit label wording in the scripts and
regenerate, since the exported glyphs are shapes rather than editable text.
`verify_vector_exports.py` checks for absent font resources, live text and
raster images, so both formats remain entirely vector graphics.

## Current plot lookup

Analysis folder names describe their contents and do not depend on figure order.
The current manuscript panel IDs and export filenames are recorded here. Paths
are relative to the repository root; `python reproduce_plots.py` rebuilds all
CLI-based panels in the table. The training folders keep their existing names.

| Current plot ID | Plot | Local source data | Reproduction entry point | Default PDF |
| --- | --- | --- | --- | --- |
| Fig2C | Normalized surface expansion | [mesh-morphometry/source_data](mesh-morphometry/source_data) | `python mesh-morphometry/reproduce_mesh_plots.py --figure 2` | [Fig2C-area_volume_evolution_control_cleanmesh_plus.pdf](mesh-morphometry/results/Fig2C-area_volume_evolution_control_cleanmesh_plus.pdf) |
| Fig4A | Model acircularity, varying thickness | [model-morphometry/ModelData/thickness_paper](model-morphometry/ModelData/thickness_paper) | [model-morphometry/reproduce_model_plots.py](model-morphometry/reproduce_model_plots.py) | [Fig4A-acircularity-varying-thickness-tol-bright.pdf](model-morphometry/results/Fig4A-acircularity-varying-thickness-tol-bright.pdf) |
| Fig4B | Experimental asphericity, BB-94 titration | [mesh-morphometry/source_data](mesh-morphometry/source_data) | `python mesh-morphometry/reproduce_mesh_plots.py --figure 4` | [Fig4B-asphericitychange_evolution_BB94_cleanmesh_plus-tol-bright.pdf](mesh-morphometry/results/Fig4B-asphericitychange_evolution_BB94_cleanmesh_plus-tol-bright.pdf) |
| Fig4E | Model bud counts | [model-morphometry/ModelData/bud_counts_summary.csv](model-morphometry/ModelData/bud_counts_summary.csv) | [model-morphometry/reproduce_model_plots.py](model-morphometry/reproduce_model_plots.py) | [Fig4E-model-bud-count-tol-bright.pdf](model-morphometry/results/Fig4E-model-bud-count-tol-bright.pdf) |
| Fig4F | Experimental bud ratio, upper panel of shared plot | [BB-94-titration/source_data](BB-94-titration/source_data) | [reproduce_combined_titration.py](reproduce_combined_titration.py), after the two analysis commands in the [titration instructions](BB-94-titration/plotting_scripts/README-titration-CLI.txt) | [Fig4FG-BB94-titration-tol-bright.pdf](BB-94-titration/results/combined/Fig4FG-BB94-titration-tol-bright.pdf) |
| Fig4G | Collagen IV intensity, lower panel of shared plot | [BB-94-titration/source_data](BB-94-titration/source_data) | Same shared renderer and preparation as Fig4F | [Fig4FG-BB94-titration-tol-bright.pdf](BB-94-titration/results/combined/Fig4FG-BB94-titration-tol-bright.pdf) |
| Fig5C | Experimental bud ratio, BB-94 washout | [BB-94-washout/source_data](BB-94-washout/source_data) | [BB-94-washout/plotting_scripts](BB-94-washout/plotting_scripts): run preparation, then either renderer; see the [washout instructions](BB-94-washout/README.txt) | [Fig5C-washout-log2-ratios-uniform-stars.pdf](BB-94-washout/results/combined/Fig5C-washout-log2-ratios-uniform-stars.pdf) |
| Fig5D | Model acircularity, full washout time course | [model-morphometry/ModelData](model-morphometry/ModelData) | [model-morphometry/reproduce_model_plots.py](model-morphometry/reproduce_model_plots.py) | [Fig5D-acircularity-full-washout.pdf](model-morphometry/results/Fig5D-acircularity-full-washout.pdf) |
| Fig5E | Model acircularity, post-washout interval | [model-morphometry/ModelData](model-morphometry/ModelData) | [model-morphometry/reproduce_model_plots.py](model-morphometry/reproduce_model_plots.py) | [Fig5E-acircularity-post-washout.pdf](model-morphometry/results/Fig5E-acircularity-post-washout.pdf) |
| Fig5F | Experimental asphericity, BB-94 washout | [mesh-morphometry/source_data](mesh-morphometry/source_data) | `python mesh-morphometry/reproduce_mesh_plots.py --figure 5` | [Fig5F-asphericitychange_evolution_washout_cleanmesh_plus.pdf](mesh-morphometry/results/Fig5F-asphericitychange_evolution_washout_cleanmesh_plus.pdf) |
| FigS1D | Segmentation benchmark: foreground-normalized error and IoU | `cellpose-SAM-training/BM-seg-benchmarking/` (downloaded from Figshare) | [cellpose-SAM-training/BM_model_benchmarking.ipynb](cellpose-SAM-training/BM_model_benchmarking.ipynb), starting at **Re-entry point: reload saved held-out masks**; details below | `FigS1D-BM_seg_model_benchmarking.pdf` via the additional export cell below |

Run `python model-morphometry/reproduce_model_plots.py` to regenerate all four
model panels together. Mesh and model outputs go
to the corresponding module's `results/` folder. Titration and washout outputs
go to `BB-94-titration/results/combined/` and `BB-94-washout/results/combined/`.
Source manifests use the current repository paths; `original_project_path`
retains the original input locations for provenance.

Manuscript plot filenames begin with the complete current panel ID. PDF, SVG
and PNG versions share a stem; the selected Tol Bright palette and marker alternatives retain their
suffixes. The combined 4F/G export starts with `Fig4FG-`. Matching numerical
reports use the same prefix. Review sheets also start with their panel IDs;
the all-panel style proofsheet is a collection overview. Original source names
remain recorded in provenance manifests.

### FigS1D: reproduce directly from the notebook

The existing notebook is the reproduction entry point; no separate CLI is
needed. Download `BM-seg-benchmarking.zip` from the `cellpose-SAM-training` item
in the [Figshare collection](https://doi.org/10.25378/janelia.c.8770956) and unzip
it inside `cellpose-SAM-training/`. The re-entry cell loads the five
`mask_gt-0.tif` through `mask_gt-4.tif` files and the five matching
`mask_model-0.tif` through `mask_model-4.tif` files. These full masks are not
included in the small repository test dataset.

Launch Jupyter from `cellpose-SAM-training/`. In
`BM_model_benchmarking.ipynb`, run cells 30-43, counting all cells from 1
(including Markdown), starting with code cell ID `a692e5e3`. This reloads the
fixed held-out split, evaluates the intended cropped regions, and regenerates
the existing vertical and horizontal plot alternatives. Stop after the two
plotting cells, IDs `7dbe859e` and `b2f9dfe1`; subsequent cells belong to image
export and model inspection. This route needs NumPy, tifffile, pandas,
Matplotlib, seaborn and Arial, without model training or inference.

The unchanged notebook's plotting cells produce these SVGs:
`jupyter_figures/BM_seg_model_benchmarking.svg` and
`jupyter_figures/BM_seg_model_benchmarking-bar.svg`. For outlined PDF and SVG
exports, plus PNG previews, run this additional cell immediately after each
plotting cell, while its `fig` and `outputFigPath` are current:

```python
from pathlib import Path
import sys
sys.path.insert(0, str(Path.cwd().parent))  # repository root
import plot_style
export_stem = Path(outputFigPath).with_suffix("")
export_stem = export_stem.with_name("FigS1D-" + export_stem.name)
plot_style.save_native_figure(fig, export_stem, dpi=600)
```

The additional cell writes `FigS1D-BM_seg_model_benchmarking` and
`FigS1D-BM_seg_model_benchmarking-bar` as PDF/SVG/PNG. It uses the existing
notebook figures and their current style, without
extracting a new CLI or changing the notebook. FigS1D remains separate from the
10 plot sets generated by the root CLI workflow. All files in the training
folders are unchanged.

## Plot sizes and style

| Panel | Full width x height (inches) | Outputs |
| --- | --- | --- |
| 2C: normalized surface expansion | 3.16 x 1.26 | [Mesh results](mesh-morphometry/results); exact frame 2.7293 x 0.96 |
| 4A: model 2D acircularity, varying z | 1.95 x 0.96 | [Model results](model-morphometry/results); Tol Bright only |
| 4B: experimental 3D asphericity | 1.95 x 0.96 | [Mesh results](mesh-morphometry/results); matched Tol Bright only |
| 4E: model bud counts | 1.32 x 0.92 | [Model results](model-morphometry/results); Tol Bright; 3 simulations per thickness |
| 4F/G: bud ratio above collagen IV intensity | 1.32 x 1.51 together | [Combined results](BB-94-titration/results/combined); Tol Bright |
| 5C: washout bud ratio | 1.00 x 1.40 | [Washout results](BB-94-washout/results/combined); two marker variants |
| 5D: full model washout time course | 1.60 x 1.24 | [Model results](model-morphometry/results) |
| 5E: model post-washout interval | 1.60 x 1.24 | [Model results](model-morphometry/results) |
| 5F: experimental washout asphericity | 3.36 x 1.24 | [Mesh results](mesh-morphometry/results) |

4A and 4B share a 110% upper y limit; 4B retains negative observations down to
-10%. Their plotting frames have the same width (40.80813036 mm, approximately
1.60662 inches) and aligned left/right edges for stacking. The final 1.5 tick
label in 4A is right-aligned inside the canvas. Their legends are raised
slightly. The 5D/5E plotting-frame height is
linked to 5F (approximately 0.9332 inches).
Figure 4B uses 0.5 pt individual curves at 30% opacity, with
opaque 1 pt means and 20%-opacity SD bands.

4E/F/G use 6 pt axis text. The horizontal `z:` and the two-line `[BB-94]` / `(µM)`
block end at the left edge of the corresponding y-axis tick labels. Intentional
left overflow is retained in SVG and shown in the review sheet, while nominal
canvas sizes remain fixed. Observations are hollow circles (2 pt in 4E/G;
1.5 pt in the denser 4F panel), with
0.35 pt outlines in the selected palette and white interiors. Black mean and
SD marks are drawn above the circles in 4E/F/G. Swarms accommodate the larger
circles; frame sizes, means, SD and endpoints are unchanged. Statistical tests
were updated to Dunnett as requested on 18 October 2026.

Figure 4E contains all 18 intended independent simulation counts: three at each
of six thicknesses, counted at tau/tau_g = 1.5. The updated table shared on
18 October 2026 completes the four previously missing runs. See the
[internal data-status note](INTERNAL_NOTES.md) and [model reproduction details](model-morphometry/README.md).

4E compares mean bud counts at each of the five thicker z values with z = 0.0005
using two-sided, single-step Dunnett tests. All groups have n = 3. Dunnett-adjusted
P values are 0.509792, 0.00141821, 0.00000391365, 0.000000841974 and
0.000000841974 (`ns`, `**`, `****`, `****`, `****`). The pooled within-group
variance is 1/3 with 12 residual degrees of freedom. Unbudded shapes count as 1.
The normal/common-variance approximation cannot be established from three runs;
this assumption is recorded explicitly in the reports.

## Selected Figure 4 colors

Paul Tol Bright is the sole Figure 4 style. Alternate black/red and rainbow
exports have been removed and are no longer generated. Output stems retain
`-tol-bright` to match existing assembly references.

4A and 4E use purple, blue, cyan, green, yellow and red for increasing thickness.
4F/G use the same order for increasing nonzero BB-94 concentration, with a gray
vehicle control. 4B maps 0, 1 and 5 uM BB-94 to the colors of z = 0.001, 0.005
and 0.01: blue, green and red. Experimental data values and panel dimensions
are unchanged by selecting this palette.

The 5C control gray (`#A9A9A9`) and BB-94 red (`#B22222`) match 5F exactly.
Control mean/SD marks are black, and 45-degree x tick labels allow a taller
frame. The canvas is 1.00 x 1.40 inches, with extra x padding beside the outer
categories. The uniform 5C version uses 1.5 pt hollow circles with 0.35 pt
colored outlines, matching 4E/F/G; mean/SD marks sit above the points. The
dataset-symbol version retains filled circles and triangles for comparison.
Palette values and provenance, including the
[Paul Tol Bright source](https://sronpersonalpages.nl/~pault/), are in `plot_palettes.json`.

The optional [review sheets](validation/review) place vector panels at their
requested physical sizes and show labels that extend beyond a panel canvas.
To rebuild them, install the optional `reportlab` dependency, then run:

```bash
python make_review_sheets.py
```

See the [current 4E/F/G review](validation/review/Fig4EFG-titration-review.pdf) and
[all-panel preview](validation/style-proofsheet.png).

## Reproduction modules

- [Mesh morphometry](mesh-morphometry/README.md): 2C, 4B and 5F from 28 tiny tables.
- [Model morphometry](model-morphometry/README.md): 4A, 4E, 5D and 5E from nine model inputs.
- [Titration instructions](BB-94-titration/plotting_scripts/README-titration-CLI.txt):
  experimental 4F/G tables and their combined renderer `reproduce_combined_titration.py`.
- [Washout instructions](BB-94-washout/README.txt): 5C preparation and both marker variants.

The two titration statistics directories retain their established CSV/report
names. The current figure exports are the combined 4F/G files; superseded
individual figure exports are not included in the main workflow.

All CLI plot paths resolve within this package, without Jupyter, manuscript
folders or network access after dependencies are installed. FigS1D instead uses
the notebook and Figshare masks described above. Reproduction replaces generated
results and preserves raw inputs. Font files are not bundled. On Windows use
`.venv\Scripts\activate` to activate the environment.

For experimental plots only on a system without TeX or Arial:

```bash
python reproduce_plots.py --skip-models --font "DejaVu Sans" --mesh-text-rendering mathtext
```

That uses approximate typography and preserves the numerical checks. Model
exports require the default Arial/Arev XeLaTeX rendering. The root requirements
are for manuscript plotting; image-processing training environments are separate.

## Model training and test data

The training folders contain code, installation instructions, and small test
crops for the two models used to process the SMG live-imaging data. The full
training datasets and trained weights are deposited in the Figshare collection
for Tauber et al., 2026 (DOI: [10.25378/janelia.c.8770956](https://doi.org/10.25378/janelia.c.8770956)).
The included test data are crops from those datasets.

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

## Bud-count and collagen-IV statistical methods

Figure 4E/F/G also have Spearman permutation trend tests reported in caption
notes rather than replacing the pairwise labels. The statistic uses all model
runs (E) or individual gland values (F/G), with ties preserved and 999,999
outcome-rank permutations (seed 20261018). Two-sided P uses absolute rho and
the nonzero `(extreme + 1)/(B + 1)` Monte Carlo correction. Holm correction
covers the three trend tests, separately from the pairwise comparison families.
Results and sample sizes are recorded in
[Fig4EFG-spearman-trends.json](validation/Fig4EFG-spearman-trends.json) and its
[readable report](validation/Fig4EFG-spearman-trends.txt).
Figure 4E/F/G use two-sided, single-step Dunnett comparisons of means against
baseline, as selected by Shaohe on 18 October 2026. E compares five thicknesses
with the thinnest z; F/G each compare six nonzero doses with DMSO (0). These are
three separate comparison families. No additional Holm adjustment is applied to
Dunnett P values, and treated groups are not compared with one another.
Horizontal bars show means and error bars show sample SD. The statistical unit
is one independent simulation (E) or gland (F/G).

Dunnett assumes independent observations and approximately normal within-group
errors with one common variance. Reports record the pooled residual variance,
degrees of freedom and, for F/G, a variance diagnostic. The root
`dunnett_statistics.py` uses deterministic conditional normal/chi-square
quadrature to evaluate the multivariate-t union tail directly, avoiding inaccurate
subtraction near one for small P values. Tests check the documented SciPy example,
the two-group pooled t limit, numerical convergence and marginal/Bonferroni bounds.

Bud counts use `log2(day1/day0)`, pairing each gland across days first. There are
110 glands from three cohorts (33, 49 and 28). The 0.005 and 0.05 uM groups have
15 glands each; other doses have 16. Two filename-marked placeholder glands are
excluded explicitly. Future-washout day2 observations are outside this endpoint.

Collagen IV uses the mean of each gland's retained line-scan AUC measurements:
716 scans, 49 glands, seven glands per dose, 12-15 scans per gland. Divide the
gland means by the mean of the seven DMSO gland means, making control mean 1.
These source tables are measurements from raw images; no denoising is applied
by the plotting workflow. The statistical unit is the gland, not the line scan.

Fig. 5 uses two-sided Welch t-tests to compare mean log2 ratios between control
and BB-94 for each of two successive intervals, with Holm correction across
those two comparisons. Welch permits unequal group variances; independent
treatment groups and approximately normal errors or an adequate mean approximation
are assumed. There are 21 glands per treatment
(8 + 13 from two datasets), with matching across time retained in each ratio.
Recorded 0/22/46 hours and annotated day0/day1/day2 are aligned by experimental
stage; ratios are not converted to rates. Symbols represent datasets when used.

Fig. 5 uses Welch-Satterthwaite degrees of freedom and retains within-gland
matching across time when calculating each ratio. No gland is duplicated within
an interval test. Holm correction allows dependence between the two interval
comparisons. Stars use adjusted p: `*` < .05; `**` < .01; `***` < .001; `****` < .0001.
All analyses treat distinct glands as independent across treatment groups.

## Plot source data and reproducibility

The original 47 experimental source tables are included byte-for-byte: 28 cleaned-mesh
measurement tables, eight Fig. 4 bud-count/key files, seven collagen IV line-scan
CSVs, and four Fig. 5 bud-count/key files. The bud-count Figure 5 ratios are rebuilt
from the raw count/key pairs; the scripts do not rely on copied intermediate
ratio exports.

Nine model inputs are also included unchanged: eight coordinate archives used
by the acircularity plots and the completed bud-count summary CSV. The CLI plotting
workflow requires no raw images, trained image-processing models, or manuscript
files. Small imaging test datasets remain in the separate training folders.

The verifier checks all 56 data-source hashes and compares regenerated gland values,
summaries and statistics to the accepted analysis. It also checks 56 mesh curves
and 7 population-SD bands, including 23 curves and three bands recolored for 4B.
The model verifier checks 11 original time courses and all 18 bud counts
against the supplied table and source-script references. Mesh time courses
use available glands at each recorded row, without interpolation or inferential
tests; see [the calculation details](mesh-morphometry/README.md#calculations-retained-from-justin).
Verification runs automatically after reproduction, or separately with:

```bash
python verify_reproduction.py
python verify_vector_exports.py --include-review  # after rebuilding review sheets
```

For unit tests:

```bash
python -m unittest discover -s BB-94-titration/plotting_scripts -p "test_*.py"
python -m unittest discover -s BB-94-washout/plotting_scripts/tests -p "test_*.py"
python -m unittest discover -s model-morphometry -p "test_*.py"
python -m unittest test_dunnett_statistics test_trend_statistics
```

Paths in generated provenance records are stored relative to this package when
running `reproduce_plots.py`. Individual CLIs may report absolute paths for their
current execution location; these are audit information, not external inputs.

See the [Figure 4 CLI instructions](BB-94-titration/plotting_scripts/README-titration-CLI.txt)
and [Figure 5 instructions](BB-94-washout/README.txt) for individual commands
and file formats.

## Method references

- [SciPy: Dunnett test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.dunnett.html)
- [SciPy: Welch t-test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html)
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

Historical presentation audits from 17 October are in `validation/history/2026-10-17/`; they describe the earlier multi-palette and incomplete-count state. Current checks are in `validation/last-verification.json`, `validation/vector-outline-checks.json`, and `validation/packaging-checks.json`.
