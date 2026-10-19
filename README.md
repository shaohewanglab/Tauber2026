# Tauber 2026 modeling

Source data and code for the Figure 2, 4 and 5 experimental and simulation
plots, the Figure S1D segmentation benchmark, and 3D-RCAN and Cellpose-SAM
model training.

## Find a plot

The links below open the manuscript PDFs. Matching SVG and PNG files are in
the same folders, alongside numerical summaries and reports.

| Panel | Plot | PDF |
| --- | --- | --- |
| Fig2C | Normalized surface expansion | [PDF](mesh-morphometry/results/Fig2C-area_volume_evolution_control_cleanmesh_plus.pdf) |
| Fig4A | Model acircularity, varying thickness | [PDF](model-morphometry/results/Fig4A-acircularity-varying-thickness-tol-bright.pdf) |
| Fig4B | Experimental asphericity, BB-94 titration | [PDF](mesh-morphometry/results/Fig4B-asphericitychange_evolution_BB94_cleanmesh_plus-tol-bright.pdf) |
| Fig4E | Model bud counts | [PDF](model-morphometry/results/Fig4E-model-bud-count-tol-bright.pdf) |
| Fig4F | Experimental bud ratio, BB-94 titration | [PDF, upper panel](BB-94-titration/results/combined/Fig4FG-BB94-titration-tol-bright.pdf) |
| Fig4G | Collagen IV intensity, BB-94 titration | [PDF, lower panel](BB-94-titration/results/combined/Fig4FG-BB94-titration-tol-bright.pdf) |
| Fig5C | Experimental bud ratio, BB-94 washout | [PDF](BB-94-washout/results/combined/Fig5C-washout-log2-ratios-uniform-stars.pdf) |
| Fig5D | Model acircularity, full washout time course | [PDF](model-morphometry/results/Fig5D-acircularity-full-washout.pdf) |
| Fig5E | Model acircularity, post-washout interval | [PDF](model-morphometry/results/Fig5E-acircularity-post-washout.pdf) |
| Fig5F | Experimental asphericity, BB-94 washout | [PDF](mesh-morphometry/results/Fig5F-asphericitychange_evolution_washout_cleanmesh_plus.pdf) |
| FigS1D | Segmentation benchmark: foreground-normalized error and IoU | [Notebook and export instructions](cellpose-SAM-training/README.md#reproduce-figs1d-segmentation-benchmark) |

Fig4F and Fig4G share one export. Fig5C also has a
[version with dataset-specific symbols](BB-94-washout/results/combined/Fig5C-washout-log2-ratios-symbols-stars-no-legend.pdf).
A [plot overview](validation/style-proofsheet.png) and
[review sheets](validation/review) are available for inspecting the exports.
PDF and SVG text is outlined, so viewing the plots requires no installed fonts.

## Reproduce the plots

Use Python 3.11 or newer. To reproduce the manuscript typography, install Arial,
TeX Live (tools `xelatex`, `xdvipdfmx`, `dvisvgm`; packages `fontspec`, `amsmath`,
`arevmath`, `pgf`, `standalone`), Poppler (`pdftoppm`) and Ghostscript (`gs`).

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce_plots.py
```

On Windows, activate the environment with `.venv\Scripts\activate`.

This regenerates the Figure 2, 4 and 5 plots, numerical summaries and statistics,
then checks source-data hashes and numerical results. All inputs are included
locally; raw images and trained models are not needed. Outputs are written to
the folders linked above. Re-running replaces generated results and preserves
source data.

FigS1D uses the [benchmarking notebook](cellpose-SAM-training/BM_model_benchmarking.ipynb)
and saved masks downloaded from Figshare. Follow the
[FigS1D instructions](cellpose-SAM-training/README.md#reproduce-figs1d-segmentation-benchmark)
to regenerate it without model training or inference.

For experimental plots on a system without TeX or Arial:

```bash
python reproduce_plots.py --skip-models --font "DejaVu Sans" --mesh-text-rendering mathtext
```

This uses approximate typography and retains the numerical checks. Model plots
require the default rendering dependencies above.

To check existing outputs without regenerating them:

```bash
python verify_reproduction.py
```

## Source data and individual analyses

Use the module instructions to reproduce selected panels or inspect their
calculations and statistical methods.

| Analysis | Panels | Source data | Instructions |
| --- | --- | --- | --- |
| Mesh morphometry | 2C, 4B, 5F | [Measurement tables](mesh-morphometry/source_data) | [Mesh analysis](mesh-morphometry/README.md) |
| Model morphometry | 4A, 4E, 5D, 5E | [Simulation inputs](model-morphometry/ModelData) | [Model analysis](model-morphometry/README.md) |
| BB-94 titration | 4F, 4G | [Bud counts and collagen IV measurements](BB-94-titration/source_data) | [Titration analysis](BB-94-titration/plotting_scripts/README-titration-CLI.txt) |
| BB-94 washout | 5C | [Bud counts](BB-94-washout/source_data) | [Washout analysis](BB-94-washout/README.txt) |
| Segmentation benchmark | S1D | Saved masks in the Figshare item `cellpose-SAM-training` | [Benchmark analysis](cellpose-SAM-training/README.md#reproduce-figs1d-segmentation-benchmark) |

Figures 4E/F/G use Dunnett comparisons against baseline; Figure 5C uses Welch
t-tests with Holm correction. Figure 4E/F/G trend tests are documented in the
[Spearman trend report](validation/Fig4EFG-spearman-trends.txt). Per-panel
statistical reports accompany the outputs in each analysis folder.

## Model training and test data

The training folders contain installation instructions, training code and small
imaging test crops. Full training datasets, benchmarking masks and trained
weights are in the [Figshare collection](https://doi.org/10.25378/janelia.c.8770956)
for Tauber et al., 2026.

| Instructions | Model and purpose | Figshare items |
| --- | --- | --- |
| [3D-RCAN training](3D-RCAN-training/README.md) | **RCAN-denoise-SMG-BM**: denoising of both live-imaging channels | `3D-RCAN-training`, `trained-models` |
| [Cellpose-SAM training](cellpose-SAM-training/README.md) | **SMG-BM-v0907**: production epithelium segmentation; **SMG-BM-15V**: benchmarking | `cellpose-SAM-training`, `trained-models` |

Each training README covers installation, a quick test and training reproduction.
Training uses separate environments; the root `requirements.txt` is for plotting.

## License

| Material | License |
| --- | --- |
| Original code, notebooks, software configuration and documentation, except the file below | [MIT](LICENSE-CODE) |
| Modified `3D-RCAN-training/apply.py` | [CC BY-NC 4.0](3D-RCAN-training/LICENSE-APPLY), noncommercial use only |
| Source and test data, derived data, plots and statistical reports | [CC BY 4.0](LICENSE-DATA) |

See [License scope and attribution](LICENSE.md) for coverage and how to credit
the data.
