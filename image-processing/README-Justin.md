# Code for analyzing 4D-microscopy data of salivary gland morphogenesis.
This folder contains scripts used for segmentation of salivary glands.

Soft Math Lab (Harvard University), in collaboration with Shaohe Wang.
Author: Justin Tauber 

## Preprocessing
- CreateClippingMask.py
- ApplyClippingMaskToBinary.py

## Automated segmentation of two-photon data into an epithelial volume and the basement membrane
- SegmentationToBinaryFile_EpitheliumBM.py

## Surface mesh generation and analysis 
- CreateMeshFromBinary.py
- AreaVolumeFromCoarseMesh_epithelial.py
- AreaVolumeFromMesh_epithelial.py

## Helper script
- ConvertBinarytoTiff.py
- ViewBinary.py
