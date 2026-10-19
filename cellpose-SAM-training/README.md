# Cellpose-SAM epithelium segmentation: training, benchmarking and test

This folder has the notebooks that trained the two Cellpose-SAM models used in the paper.
- **SMG-BM-v0907**: the production model, used for all `*_seg_model0907` segmentations.
- **SMG-BM-15V**: the benchmarking model. It was trained on 15 of the 20 annotated volumes
  and scored on the other 5.

The folder also has a small test volume and a test script. The full training data, the
benchmarking masks and the trained weights are in the Figshare collection
(DOI: [10.25378/janelia.c.8770956](https://doi.org/10.25378/janelia.c.8770956)):

| Figshare item | Contents |
|---|---|
| `cellpose-SAM-training` | Training volumes, manual slices, benchmarking masks, test volumes and these notebooks |
| `trained-models` | `cellpose-SMG-BM-v0907.zip`, `cellpose-SMG-BM-15V.zip` (trained weights) |

## Files

| File | Contents |
|---|---|
| `basement_membrane_training.ipynb` | Trains SMG-BM-v0907 on all 20 annotated volumes plus manually annotated low-SNR slices, then applies it in 3D |
| `BM_model_benchmarking.ipynb` | Trains SMG-BM-15V on 15 volumes and scores it on the 5 held-out volumes (voxel error rate and IoU). It can be re-run from the saved masks |
| `test_cellpose.py` | Smoke test on `test_data/` (see below) |
| `test_data/` | A 147 x 160 x 160 (z, y, x) crop of a volume held out from SMG-BM-15V training (`5_SMG4_t090`, dataset 241104-SMG4): `test-c1-rcan.tif` (Col4a2-mNG), `test-c2-rcan.tif` (PNA-Cy3), RCAN-denoised 16-bit; `test-mask-annotated.tif`, annotation; `crop.json`, crop position. Voxel size 0.8546 x 0.8546 x 1.0 um |
| `cellpose_installation.md`, `cellpose_requirements.txt` | Installation: the Cellpose instructions, or the exact environment we used (pinned in `cellpose_requirements.txt`) |

## Install

Install Cellpose following the Cellpose repository
(https://github.com/MouseLand/cellpose#installation), plus a few extra packages. The exact
environment we used was Cellpose 4.0.4, Python 3.12 and PyTorch 2.7.1 with CUDA 12.6. Both
options are described in `cellpose_installation.md`.

## Quick test

From this folder:

```bash
# 1. inference with the trained model: unzip cellpose-SMG-BM-v0907.zip from the Figshare
#    item trained-models and put the model file at models/SMG-BM-v0907
python test_cellpose.py --model models/SMG-BM-v0907

# 2. training check: fine-tunes a throw-away model for 5 epochs on slices of the test volume
python test_cellpose.py --train
```

The script segments the test volume in 3D with the settings used in the paper. It
thresholds the cell probability at 0, fills holes and keeps the largest object. It then
prints the voxel error rate and IoU against the annotation and writes
`test_data/test-seg.tif`. Add `--cpu` to run without a GPU (slow).

Expected result with SMG-BM-v0907: a voxel error rate of about 0.044 and an IoU of about
0.956 (Cellpose 4.0.9, CPU). Small differences between Cellpose versions and between GPU
and CPU are normal. The `--train` check only shows that training runs; its scores are not
meaningful.

With `--train`, a few "Mean of empty slice" warnings during training are harmless; they
come from slices without annotated tissue. A PyTorch "Sparse invariant checks are
implicitly disabled" warning during segmentation is also harmless. On our GPU the `--train`
check gave an error rate of about 0.16 and an IoU of about 0.86 after 5 epochs. After training, the script loads the saved
model (about 1.2 GB, written to `models/test-finetune`) and segments the volume in 3D.
Cellpose shows no progress during that step. Delete `models/` and `test_data/test-seg.tif`
after testing.

## Reproducing the training

1. Download the Figshare item `cellpose-SAM-training` and unzip its zips in this folder.
   This gives `SMG-BM-seg-training/`, `train_manual/`, `BM-seg-benchmarking/` and
   `BM-img-to-test/` next to the notebooks.
2. Start Jupyter here and run the notebooks. All data paths are relative to this folder.
   Trained models are written to `models/<model_name>` and figures to `jupyter_figures/`.

Training samples random slices, and the 15/5 benchmarking split is not seeded; the split
used for the saved masks is recorded in the notebook. Re-training therefore gives a
similar but not identical model.

## Reproduce FigS1D: segmentation benchmark

The existing notebook is the reproduction entry point; no separate CLI is
needed. Download `BM-seg-benchmarking.zip` from the `cellpose-SAM-training` item
in the [Figshare collection](https://doi.org/10.25378/janelia.c.8770956) and unzip
it inside this folder. The re-entry cell loads the five
`mask_gt-0.tif` through `mask_gt-4.tif` files and the five matching
`mask_model-0.tif` through `mask_model-4.tif` files. These full masks are not
included in the small repository test dataset.

Launch Jupyter from this folder. In
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
10 plot sets generated by the root CLI workflow.

## Licenses

The notebooks and test script are covered by the repository's
[MIT License](../LICENSE-CODE). The bundled `test_data/`, including crop
metadata, are licensed under [CC BY 4.0](../LICENSE-DATA). Cellpose and other
third-party dependencies retain their own licenses; see
[repository licensing](../LICENSE.md).
