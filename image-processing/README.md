# Image processing, clipping and meshing

Reference scripts shared by Justin Tauber, Soft Math Lab (Harvard University),
in collaboration with Shaohe Wang. The eight Python files are preserved as
supplied, as a starting point for adapting the rule-based segmentation and
clipping/meshing tools to other data. Edit their dataset-specific settings
before use. The supplied overview is retained as [README-Justin.md](README-Justin.md),
and [provenance.json](provenance.json) records the original file hashes.

## Scripts

| Script | Purpose |
| --- | --- |
| `SegmentationToBinaryFile_EpitheliumBM.py` | Rule-based epithelial and basement-membrane segmentation: Hessian enhancement, morphological reconstruction, watershed refinement and pruning; uses manual seeds and optionally the previous frame |
| `CreateClippingMask.py` | Draw polygons in napari to define regions to exclude, such as the primary duct; saves a 3D `reference_mask.tif` |
| `ApplyClippingMaskToBinary.py` | Remove the clipping-mask region from Cellpose segmentation TIFFs, retain the largest component and save binary `.npz` volumes |
| `CreateMeshFromBinary.py` | Marching cubes, largest-component selection, Taubin smoothing and decimation; saves coarse `.ply` meshes |
| `AreaVolumeFromCoarseMesh_epithelial.py` | Measure coarse meshes in `Results/` |
| `AreaVolumeFromMesh_epithelial.py` | Measure cleaned meshes in `Results_MeshVoxelization/` |
| `ConvertBinarytoTiff.py` | Export binary `.npz` volumes as calibrated ImageJ TIFFs |
| `ViewBinary.py` | Inspect binary volumes in napari |

Clipping acts on the segmented volume **before meshing**. These scripts do not
implement interactive cutting of an existing mesh.

## Setup and adaptation

Use a separate Python environment from the plotting and model-training
environments. From this folder:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate with `.venv\Scripts\activate`. The requirements list the
imported packages, with `napari[all]` supplying a Qt backend as described in the
[napari installation guide](https://napari.org/stable/getting_started/installation.html).
Versions are not pinned: Justin's original environment was not supplied, and
this deposit has been checked for source integrity and Python syntax, without
an end-to-end run.

Run the files as standalone scripts. Edit the settings near the main block
(or near the top of the two measurement scripts): paths, dataset names,
filename prefixes, frame ranges and voxel calibration. Some paths still use
Justin's `/mnt/d/...` working directory. The measurement scripts assume five
minutes per frame; change this for other acquisitions.

Relative paths are resolved from the working directory. Create the parent
`Data/`, `Annotation/` and `Results/` directories there before running, and
invoke a script by its full path if working outside this folder. Local data
directories within this folder are ignored by Git.

Arrays are in `(z, y, x)` order. In the supplied meshing and TIFF-export scripts,
the variable named `voxel_size_x` actually supplies the stack's z spacing;
check how calibration is used rather than relying on variable names alone.

## Rule-based segmentation

This is the draft-segmentation approach used before manual curation of the
reference annotations. The production image series were segmented with
[Cellpose-SAM](../cellpose-SAM-training/README.md). Denoising is described in
[3D-RCAN training and application](../3D-RCAN-training/README.md).

`SegmentationToBinaryFile_EpitheliumBM.py` expects paired, denoised, 16-bit TIFF
stacks under `Data/<dataset>/`, with filenames set in the script. Prepare your
own annotations under `Annotation/<SaveFolder>/`:

- `reference_mask.tif`: clipping/exclusion mask, matching the stack dimensions.
- `rough_mask_frame_000.tif`: a rough foreground region used to constrain the
  initial segmentation.
- `interior_points_frame_<frame>.csv`: interior seed coordinates, with a header
  and columns `z,y,x`. Supply seeds for the first processed frame. The supplied
  loader expects multiple point rows; later seed files replace the previous
  coordinates when present.

The script saves `thresholded_epithelial_<frame>.npz` and
`thresholded_membrane_<frame>.npz`, each with an array named `data`.
The rough-mask and seed-point annotations are user inputs; tools for creating
them are not included in this folder.

## Clipping and meshing Cellpose outputs

After editing the settings in each script, run them in this order:

1. `CreateClippingMask.py`: displays the first and last segmentation projections.
   Draw the regions to remove in the `polygon_regions` layer, then close napari
   to save the polygons extruded through z as `reference_mask.tif`.
2. `ApplyClippingMaskToBinary.py`: applies that exclusion mask to the configured
   segmentation TIFFs and writes `Results/<SaveFolder>/thresholded_epithelial_<frame>.npz`.
3. `CreateMeshFromBinary.py`: converts those volumes to
   `Results/<folder>/mesh_epithelial_<frame>.ply`.
4. `AreaVolumeFromCoarseMesh_epithelial.py`: measures the coarse meshes.

`ConvertBinarytoTiff.py` and `ViewBinary.py` provide export and inspection of
the binary volumes. Coordinate the dataset settings across scripts; their
supplied defaults refer to different datasets.

## Relationship to PaperData and Figshare

The full `PaperData` folder will be deposited on Figshare. The paper's
[Figshare collection](https://doi.org/10.25378/janelia.c.8770956) already provides
training data and models; a direct link to the `PaperData` item can be added
when available. Dataset-specific clipping masks, manual annotations, images
and meshes are kept with the data deposit, outside this code folder.

| PaperData location | Relationship to these scripts |
| --- | --- |
| `overview_data_paper.xlsx` | Dataset index and selection notes |
| `Annotation/<dataset>/` | Clipping masks made with `CreateClippingMask.py` |
| `Results/<dataset>_cp_model0907/` | Clipped epithelial `.npz` volumes, coarse `.ply` meshes and `volume_area_*.txt` tables |
| `Results_MeshVoxelization/<dataset>_cp_model0907_v2/` | Cleaned meshes and tables measured with `AreaVolumeFromMesh_epithelial.py` |
| `Results_SphericalHarmonics/`, plus `spherical_harmonics/` subfolders in some cleaned-mesh datasets | Downstream `.mat` coefficients and `.stl` reconstructions; not generated by these eight scripts |

The supplied scripts cover coarse mesh generation and measurement of already
cleaned meshes. The intervening voxelization/repair workflow is not included.
The rule-based script's rough masks and seed CSVs are also absent from the
shared `PaperData` folder; it primarily holds the production Cellpose results.

Tables contain `frame`, `time_min`, `volume_um3` and `surface_area_um2`.
Non-watertight meshes receive a `nan` volume. The 28 cleaned-mesh tables used
for Figures 2C, 4B and 5F are already included in
[mesh-morphometry/source_data](../mesh-morphometry/source_data), with their
original `PaperData` paths and hashes in its manifest. See the
[mesh-morphometry instructions](../mesh-morphometry/README.md) to reproduce
those plots directly from the tables.

## Attribution and license

Justin's author headers and references to RACE and spheresDT/Mpacts-PiCS are
retained in the source. The scripts and documentation are covered by the
repository's [MIT License](../LICENSE-CODE); dependencies retain their own
licenses. See [repository licensing](../LICENSE.md).
