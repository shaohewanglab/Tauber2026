# Cellpose-SAM installation instructions

## Recommended: follow the Cellpose instructions

Install Cellpose following the instructions in the Cellpose repository:
https://github.com/MouseLand/cellpose#installation

Those instructions are kept up to date for current Python, PyTorch and CUDA versions and
for each operating system. Cellpose-SAM needs an NVIDIA GPU for practical speed; on the CPU
it is very slow. The notebooks and `test_cellpose.py` work with Cellpose 4 (Cellpose-SAM).

Then add the few extra packages used by the notebooks and the test script:

```bash
pip install jupyter matplotlib natsort opencv-python-headless pandas scipy seaborn tifffile
```

## The environment we used

The notebooks were run with Cellpose 4.0.4, Python 3.12.11 and PyTorch 2.7.1 with CUDA 12.6,
on Windows. To reproduce it:

```bash
conda create -y -n cellpose python=3.12
conda activate cellpose

# PyTorch with CUDA 12.6 (see https://pytorch.org for other CUDA versions or platforms)
pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu126

# Cellpose 4.0.4 and the other packages used by the notebooks and test_cellpose.py
pip install -r cellpose_requirements.txt
```

Check that PyTorch sees the GPU:

```bash
python -c "import torch; print(torch.cuda.is_available())"
```
