"""
Script for creating a clipping mask using Napari.
Used for e.g. clipping the primary duct.
Soft Math Lab (Harvard University)
Author: Justin Tauber

Collaboration with Shaohe Wang
"""

import time

import numpy as np
from pathlib import Path

import tifffile as tiff

import napari

if __name__ == "__main__":

	FileFolder = "250916-SMG2-BB94-washout"
	FileName = "250916-Col4a2-mNG-BB94-washout-003_2_t"
	ImageName = FileName
	FilePath1 = f"./Data/{FileFolder}/{FileName}split_c1_rcan/"
	FilePath2 = f"./Data/{FileFolder}/{FileName}split_c2_rcan/"
	FilePath3 = f"./Data/{FileFolder}/{FileName}split_seg_model0907/"

	AnnotationPath = f"./Annotation/{FileFolder}/"
	annotationDir = Path(AnnotationPath)
	annotationDir.mkdir(exist_ok=True)

	# Range for analysis			
	start = 0
	stop = 218

##########################################################
# LOAD-FRAME
##########################################################
global_start_time = time.time()

# Get the last image and compute z projection
i = stop
print(f"\nLoading frame {i}/{stop}")

tif_data3 = tiff.imread(f"{FilePath3}{ImageName}{i:03d}-seg.tif")					
image3 = tif_data3[:,:,:] / (2**16 - 1)
z_mean_image_stop = image3.mean(axis=0)

# Get the first image and compute z projection
i = start
print(f"\nLoading frame {i}/{stop}")
tif_data3 = tiff.imread(f"{FilePath3}{ImageName}{i:03d}-seg.tif")					
image3 = tif_data3[:,:,:] / (2**16 - 1)
z_mean_image_start = image3.mean(axis=0)

zdim = np.shape(image3)[0]
xdim = np.shape(image3)[2]
ydim = np.shape(image3)[1]

# Create an empty viewer
viewer = napari.Viewer()

# Show 2D image
viewer.add_image(z_mean_image_start, name='start', colormap='magenta', blending='additive')
viewer.add_image(z_mean_image_stop, name='stop', colormap='green', blending='additive')

# Add shapes layer for drawing polygons
viewer.add_shapes(name='polygon_regions', shape_type='polygon')

napari.run()

if 'polygon_regions' in viewer.layers:

	layer = viewer.layers['polygon_regions']
	if len(layer.data) > 0:

		# load masks and combine into 1 image
		masks = layer.to_masks(mask_shape=(ydim, xdim))
		mask_2d = masks.any(axis=0)

		# Convert into 3D image of same dimensions as input
		mask_3d = np.repeat(mask_2d[None, :, :], zdim, axis=0)
		tiff.imwrite(annotationDir / "reference_mask.tif", mask_3d.astype(np.uint8))