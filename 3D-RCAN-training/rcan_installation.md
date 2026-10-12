# 3D-RCAN installation instructions

These steps work on Windows and Linux. They have not been tested on macOS, where it is
probably not worth trying because there is no CUDA support.

This is an expanded version of the installation instructions in the
[3D-RCAN repository](https://github.com/AiviaCommunity/3D-RCAN).

## 1. Get the 3D-RCAN code

3D-RCAN itself needs no installation; its scripts are run from the repository folder.

```bash
git clone https://github.com/AiviaCommunity/3D-RCAN.git
```

(or download the repository as a zip from GitHub and unzip it).

## 2. Create the conda environment

```bash
# Python 3.7 is needed because tensorflow_gpu 1.13.1 is only available for older Python versions
conda create -y -n rcan python=3.7

conda activate rcan
conda install -y cudatoolkit=10.0 cudnn=7.6.5

# install the requirements listed in rcan_requirements.txt (from this Figshare item)
pip install -r rcan_requirements.txt

# downgrade protobuf for compatibility with tensorflow_gpu 1.13.1
conda install -y protobuf=3.20.1
```

Contents of `rcan_requirements.txt`:

```text
colorama; sys_platform=='win32'
h5py<3.0
jsonschema
keras==2.2.4
numexpr
numpy>=1.18,!=1.19.4
tensorflow_gpu==1.13.1
tifffile
tqdm>=4.41.0
```

## 3. Use the modified apply.py

Copy `apply.py` and `voxel_size_utils.py` from this Figshare item into the 3D-RCAN
repository folder, replacing its `apply.py`. The two files must sit next to each other
and next to the `rcan` package folder. The modified `apply.py` keeps the voxel size and
unit of the input in the output TIFF and by default writes only the restored image.

## 4. Run

With the `rcan` environment active, call the scripts by their path in the 3D-RCAN folder.
Python finds the `rcan` package next to the script, so you can run from any folder. Data
paths, including those in the config file, are relative to the folder you run from.

```bash
# training: run from the unzipped Figshare item folder, where config_Col4a2-mNG.json and
# RCAN-training-data-SMG/ are (see README.md)
python <3D-RCAN folder>/train.py -c config_Col4a2-mNG.json -o RCAN-denoise-SMG-BM

# denoising a folder of TIFF stacks, 16-bit output
python <3D-RCAN folder>/apply.py -m RCAN-denoise-SMG-BM -i <input folder> -o <output folder> -b 16
```
