# 3D-RCAN denoising: training, application and test

This folder has the scripts and configuration for the 3D-RCAN model RCAN-denoise-SMG-BM. It
was trained on paired low- and high-SNR 3D images of Col4a2-mNG in fixed SMG explants and
used to denoise both channels of the live-imaging data. A small test set is included. The
full training data and the trained weights are in the Figshare collection
(DOI: [10.25378/janelia.c.8770956](https://doi.org/10.25378/janelia.c.8770956)):

| Figshare item | Contents |
|---|---|
| `3D-RCAN-training` | `RCAN-training-data-SMG.zip` (expanded set: 34 training and 6 validation pairs; manuscript model used the 20 + 6 September pairs) and the files in this folder |
| `trained-models` | `RCAN-denoise-SMG-BM.zip` (weights, training config, TensorBoard log) |

## Files

| File | Contents |
|---|---|
| `apply.py` | Minimally modified 3D-RCAN `apply.py`. It keeps the ImageJ voxel size and unit of the input in the output, and by default writes only the restored image (`--output_raw` and `--output_gt` add the input and ground truth as channels). The changes are listed in its header |
| `voxel_size_utils.py` | Reads and writes the ImageJ voxel size; imported by `apply.py` |
| `config_Col4a2-mNG.json` | Training configuration of RCAN-denoise-SMG-BM (data paths relative to the unzipped Figshare data) |
| `config_test.json` | Same model settings, scaled down for a quick training check on `test_data/` (4 epochs x 5 steps, 32 x 64 x 64 patches; at least 4 epochs are needed because 3D-RCAN halves the learning rate every epochs / 4 epochs; `intensity_threshold` 0, because the small crops have too little bright signal for 3D-RCAN's default patch check) |
| `test_data/` | 3 pairs (2 training, 1 validation) cropped to 61 x 128 x 128 (z, y, x) from the training data: `*_raw/*_c1.tif` low-SNR input, `*_gt/*_c2.tif` high-SNR target; 16-bit, 0.752 x 0.752 x 1.0 um, ImageJ-calibrated |
| `rcan_installation.md`, `rcan_requirements.txt` | Installation |

## Install

See `rcan_installation.md`. It covers cloning 3D-RCAN, creating the conda environment, and
copying `apply.py` and `voxel_size_utils.py` into the 3D-RCAN folder in place of its
`apply.py`.

## Quick test

From this folder, with the `rcan` environment active. Run the two commands one after the
other, not at the same time on one GPU:

```bash
# 1. training check (a few minutes on a GPU; the model is not useful)
python <3D-RCAN folder>/train.py -c config_test.json -o test_model

# 2. denoising with the trained model (unzip RCAN-denoise-SMG-BM.zip from the Figshare item
#    trained-models here); writes raw, restored and ground truth as three channels
python <3D-RCAN folder>/apply.py -m RCAN-denoise-SMG-BM -i test_data/validate_raw -g test_data/validate_gt -o test_output -b 16 --output_raw --output_gt
```

The output TIFF in `test_output/` should carry the input's voxel size (0.752 x 0.752 x 1.0
um). The deposited `*_rcan` stacks were made with `-b 16` and the default restored-only
output.

## Reproducing the training

Download `RCAN-training-data-SMG.zip` from the Figshare item `3D-RCAN-training`. For the
manuscript training setup, extract only the TIFF files with names beginning `240906`
into this folder, preserving the archive's directory structure. This gives 20 training
pairs and 6 validation pairs in the directories specified by `config_Col4a2-mNG.json`.
The 14 additional training pairs beginning `241111` were added later and were not used
for the manuscript checkpoint. Then run:

```bash
python <3D-RCAN folder>/train.py -c config_Col4a2-mNG.json -o RCAN-denoise-SMG-BM
```

The configuration specifies 300 epochs x 400 steps and a patch size of 64 x 128 x 128
(z, y, x). For the manuscript model, 3D-RCAN automatically reduced the patch size to
61 x 128 x 128 to fit the original training and validation volumes. The run completed
192 epochs; the deposited weights are from epoch 165, which had the lowest validation
loss and was used to process the manuscript data.

Using the full expanded archive (34 training and 6 validation pairs) instead reduces
the patch size to 51 x 128 x 128 because the additional November volumes have 51 slices.
This describes retraining on the expanded dataset, not the deposited manuscript model.

3D-RCAN prints "Failed to sample a valid patch" from time to time during this training. This
is expected with these sparse basement-membrane images. A training patch is accepted only if
at least half of its target voxels are above 0.25 after normalization, which very few
patches meet. After 512 failed tries the last random patch is used anyway, as in the
original training. The test config turns the check off (`intensity_threshold` 0), so the
quick test should run without these warnings.

## Licenses

`apply.py` is adapted from 3D-RCAN (https://github.com/AiviaCommunity/3D-RCAN, Copyright
2021 SVision Technologies LLC). Like 3D-RCAN, it is licensed under CC BY-NC 4.0
(https://creativecommons.org/licenses/by-nc/4.0/), which does not allow commercial use.
`voxel_size_utils.py` is new code under the MIT License (text in the file).
See [LICENSE-APPLY](LICENSE-APPLY) for the modified script's license notice and
[repository licensing](../LICENSE.md) for the overall scope. The bundled
`test_data/` are licensed under [CC BY 4.0](../LICENSE-DATA).
