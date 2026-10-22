"""
Algorithm for applying a clipping mask to 4D microscopy data
Soft Math Lab (Harvard University)
Author: Justin Tauber

Collaboration with Shaohe Wang

"""

import time
import numpy as np
from pathlib import Path

import tifffile as tiff
from skimage import measure

##################################################################
# Functions for segmentation
##################################################################

def largest_component(image):
	# Label connected components
	labels = measure.label(image, connectivity=3)

	# Find the largest component
	props = measure.regionprops(labels)
	largest_component = max(props, key=lambda x: x.area)

	# Create a mask for the largest component
	largest_component_mask = labels == largest_component.label

	return largest_component_mask

##################################################################
# Main loop
##################################################################

if __name__ == "__main__":
	FileFolder = "250916-SMG8-BB94-washout"
	SaveFolder = "250916-SMG8-BB94-washout_cp_model0907"
	FileName = "250916-Col4a2-mNG-BB94-washout-003_8_t"
	ImageName = FileName
	FilePath1 = f"./Data/{FileFolder}/{FileName}split_c1_rcan/"
	FilePath2 = f"./Data/{FileFolder}/{FileName}split_c2_rcan/"
	CellPosePath = f"./Data/{FileFolder}/{FileName}split_seg_model0907/" #

	ResultPath = f"./Results/{SaveFolder}/"
	ResultDir = Path(ResultPath)
	ResultDir.mkdir(exist_ok=True)

	AnnotationPath = f"./Annotation/{FileFolder}/"

	# Dimensions of the voxels
	voxel_size_x = 1.0
	voxel_size_y = 1.3811
	voxel_size_z = 1.3811

	# Spacing used for plotting in Napari
	original_spacing = np.array([voxel_size_x, voxel_size_y, voxel_size_z])
	normalized_spacing = original_spacing / 1.0

	# Range for analysis			
	start = 0
	stop = 218
	count = 0
	step = 1

	reference_path = f"{AnnotationPath}/reference_mask.tif"
	reference_mask = tiff.imread(reference_path).astype(bool)

	for i in range(start,stop+1,step):

		##########################################################
		# LOAD-FRAME
		##########################################################

		print(f"\nLoading frame {i}/{stop} ({count})")
		global_start_time = time.time()

		##########################################################
		# LOAD ANNOTATION AND MASKS
		##########################################################
		try:
			gland = tiff.imread(f"{CellPosePath}{ImageName}{i:03d}-seg.tif").astype(bool)
		except FileNotFoundError:
			print(f"Could not find frame {i}")
			continue
		gland_combined = np.copy(gland)
		gland_combined[reference_mask] = False

		gland_combined = largest_component(gland_combined)

		##########################################################
		# SAVE RESULTS
		##########################################################

		np.savez_compressed(f"./Results/{SaveFolder}/thresholded_epithelial_{i:03d}", data=gland_combined)

		global_end_time = time.time()
		elapsed_time = global_end_time - global_start_time
		print(f"frame time:\t{elapsed_time:.2f}")

		count +=1