"""
Segmentation algorithm for Salivary glands
Soft Math Lab (Harvard University)
Author: Justin Tauber

Collaboration with Shaohe Wang

The segmentation implemented here is inspired by two segmentation methods:
Real-time Accurate Cell-shape Extractor (RACE) https://doi.org/10.1016/j.devcel.2015.12.028
spheresDT/Mpacts-PiCS:(https://academic.oup.com/bioinformatics/article/37/24/4851/6331998)

# The algorithm:
1. reads in tiff stacks, assuming they carry denoised data (e.g. using RCAN)
2. Preprocessing: rescaling image intensity
3. Filtering: Hessian objectness filter for enhancing platelike structures
4. Segment membrane based on xz and yz slices.
5. Segment salivary gland based on closed membrane segmentation
6. Optimize gland segmentation using watershed. 
   Using initial segmentation and boundary as seed, excluding fine membrane segmentation.
7. Postprocessing: remove thin components

Segmentation of previous frames is used in guiding segmentation in the next time point.

"""

from joblib import Parallel, delayed
import time

import numpy as np
from scipy import ndimage as ndi

from pathlib import Path

import tifffile as tiff

from skimage.filters import threshold_otsu
from skimage import morphology
from skimage import measure
from skimage import filters
from skimage import segmentation

import itk
import napari

###############################################
# Functions for visualization
###############################################

def showNapari(image, spacing):
	aspect_xz = spacing[0]/spacing[1]
	aspect_yz = spacing[0]/spacing[2]

	# Create an empty viewer
	viewer = napari.Viewer()
	
	# Add a new layer containing an image
	viewer.add_image(image)
	
	# Remove all layers to start from scratch
	for l in viewer.layers:
		viewer.layers.remove(l)
		
	# add the image again with scaling
	viewer.add_image(image, scale=(aspect_xz, 1, 1))

	napari.run()

def showNapari2D(image):

	# Create an empty viewer
	viewer = napari.Viewer()
	
	# Add a new layer containing an image
	viewer.add_image(image)
	
	napari.run()

def showNapari2Im(image1, image2, spacing):
	aspect_xz = spacing[0] / spacing[1]
	aspect_yz = spacing[0] / spacing[2]

	# Create an empty viewer
	viewer = napari.Viewer()

	# Add the first image layer (using the red channel)
	viewer.add_image(image1, scale=(aspect_xz, 1, 1), colormap='magenta', blending='additive', name='Channel 1')

	# Add the second image layer (using the blue channel)
	viewer.add_image(image2, scale=(aspect_xz, 1, 1), colormap='green', blending='additive', name='Channel 2')

	napari.run()

def showNapari2ImPoints(image1, image2, spacing, interior_points, exterior_points):
	aspect_xz = spacing[0] / spacing[1]
	aspect_yz = spacing[0] / spacing[2]

	viewer = napari.Viewer()

	viewer.add_image(image1, scale=(aspect_xz, 1, 1), colormap='magenta', blending='additive', name='Channel 1')
	viewer.add_image(image2, scale=(aspect_xz, 1, 1), colormap='green', blending='additive', name='Channel 2')

	if len(interior_points) > 0:
		scaled_interior = interior_points * np.array([aspect_xz, 1, 1])
		viewer.add_points(scaled_interior, name='interior_points', face_color='blue', size=4)

	if len(exterior_points) > 0:
		scaled_exterior = exterior_points * np.array([aspect_xz, 1, 1])
		viewer.add_points(scaled_exterior, name='exterior_points', face_color='yellow', size=4)
	napari.run()

def showNapari3Im(image1, image2, image3, spacing):
	aspect_xz = spacing[0] / spacing[1]
	aspect_yz = spacing[0] / spacing[2]

	# Create an empty viewer
	viewer = napari.Viewer()

	# Add the first image layer (using the red channel)
	viewer.add_image(image1, scale=(aspect_xz, 1, 1), colormap='red', blending='additive', name='Channel 1')

	# Add the second image layer (using the blue channel)
	viewer.add_image(image2, scale=(aspect_xz, 1, 1), colormap='blue', blending='additive', name='Channel 2')

	# Add the second image layer (using the blue channel)
	viewer.add_image(image3, scale=(aspect_xz, 1, 1), colormap='green', blending='additive', name='Channel 3')


	napari.run()

def showNapari2D2Im(image1, image2, spacing):

	aspect_xz = spacing[0] / spacing[1]

	# Create an empty viewer
	viewer = napari.Viewer()
	
	# Add a new layer containing an image
	viewer.add_image(image1, scale=(aspect_xz,1), colormap='red', blending='additive', name='Channel 1')
	viewer.add_image(image2, scale=(aspect_xz,1), colormap='blue', blending='additive', name='Channel 2')	
	napari.run()

###############################################
# Functions for segmentation
###############################################

def vesselness_downhill_binarization(image, global_threshold=None, membrane_acceptance_level=100):

	# Downhill filtering parameters
	if global_threshold == None:
		l_otsu = threshold_otsu(image)
	else:
		l_otsu = global_threshold

	l_abs_tresh_slice_seed = l_otsu
	l_abs_thresh_slice = l_abs_tresh_slice_seed / membrane_acceptance_level

	# Select the seed
	seed = np.copy(image)
	seed[image < 2*l_otsu] = 0

	# Apply morphological reconstruction
	mask = np.where(image > l_abs_thresh_slice, image, 0)
	a_threshold = morphology.reconstruction(seed=seed, mask=mask)
	a_threshold = a_threshold >= l_abs_thresh_slice

	return a_threshold

def plateness_computation_spacing(image, spacing):
	"""
	Apply the Hessian-based objectness filter for membrane enhancement,
	adjusted for non-cubic voxels with specified spacing.

	Parameters:
		image (numpy.ndarray): Input 3D grayscale image (scaled 0 to 1).
		spacing (tuple): Tuple specifying voxel spacing (z, y, x).
		
	Returns:
		numpy.ndarray: Output image with plateness enhancement.
	"""
	Dimension = 3
	PixelType = itk.D
	ImageType = itk.Image[PixelType, Dimension]
	HessianPixelType = itk.SymmetricSecondRankTensor[itk.D, Dimension]
	HessianImageType = itk.Image[HessianPixelType, Dimension]

	# Convert NumPy array to ITK image
	itk_image = itk.GetImageFromArray(image.astype(np.float64))
	
	# Set the voxel spacing to account for non-cubic voxels
	itk_image.SetSpacing(spacing)

	# Set up the objectness filter
	objectness_filter = itk.HessianToObjectnessMeasureImageFilter[HessianImageType, ImageType].New()
	objectness_filter.SetBrightObject(True)
	objectness_filter.SetScaleObjectnessMeasure(False)
	objectness_filter.SetAlpha(1.0)
	objectness_filter.SetBeta(1.0)
	objectness_filter.SetGamma(0.05) # 0.01 0.05
	objectness_filter.SetObjectDimension(2)  # For plane-like structures

	# Set the multiscale Hessian filter
	multi_scale_filter = itk.MultiScaleHessianBasedMeasureImageFilter[ImageType, HessianImageType, ImageType].New()
	multi_scale_filter.SetInput(itk_image)
	multi_scale_filter.SetHessianToMeasureFilter(objectness_filter)
	multi_scale_filter.SetSigmaMinimum(0.5)
	multi_scale_filter.SetSigmaMaximum(1.5)
	multi_scale_filter.SetNumberOfSigmaSteps(3)

	# Rescale output intensity to the range of unsigned short (16-bit)
	OutputPixelType = itk.US
	OutputImageType = itk.Image[OutputPixelType, Dimension]

	rescale_filter = itk.RescaleIntensityImageFilter[ImageType, OutputImageType].New()
	rescale_filter.SetInput(multi_scale_filter)

	# Get the result
	rescale_filter.Update()
	result_image = rescale_filter.GetOutput()
	plateness = itk.GetArrayFromImage(result_image)

	return plateness

def largest_component(image):
	# Label connected components
	labels = measure.label(image, connectivity=3)

	# Find the largest component
	props = measure.regionprops(labels)
	largest_component = max(props, key=lambda x: x.area)

	# Create a mask for the largest component
	largest_component_mask = labels == largest_component.label

	return largest_component_mask

def overlapping_components(image, mask):
	# Label connected components
	labels = measure.label(image, connectivity=3)

	touching_labels = np.unique(labels[image & mask])
	overlapping_mask = np.isin(labels, touching_labels)

	return overlapping_mask

def overlapping_components_with_points(image, points):
	# image: binary image (e.g. candidate regions)
	# points: (n, 3) array of [z, y, x] interior points
	shape = np.shape(image)

	# Convert points to mask
	point_mask = np.zeros(shape, dtype=bool)
	for z, y, x in points.astype(int):
		if 0 <= z < shape[0] and 0 <= y < shape[1] and 0 <= x < shape[2]:
			point_mask[z, y, x] = True

	# Label components
	labels = measure.label(image, connectivity=3)

	# Find labels that touch the points
	touching_labels = np.unique(labels[point_mask])
	touching_labels = touching_labels[touching_labels != 0]  # Exclude background
	overlapping_mask = np.isin(labels, touching_labels)

	return overlapping_mask

def largest_component2D(image):
	# Label connected components
	labels = measure.label(image)

	# Find the largest component
	props = measure.regionprops(labels)
	largest_component = max(props, key=lambda x: x.area)

	# Create a mask for the largest component
	largest_component_mask = labels == largest_component.label

	return largest_component_mask

def freefloating_components(image):
	# Label connected components of zero value
	labels = measure.label(image, connectivity=1)

	# Get the shape of the image
	rows, cols = image.shape

	# Find components touching the boundary
	boundary_labels = np.unique(
		np.concatenate([
			labels[0, :], labels[rows - 1, :], labels[:, 0], labels[:, cols - 1]
		])
	)

	# Create a mask for components not touching the boundary
	mask = np.isin(labels, boundary_labels, invert=True)	
	return image & mask

def enclosed_black_regions(red_mask, blue_mask):
	# Step 1: Identify black regions in the red mask (gland)
	inverted_red = ~red_mask
	
	# Label connected components in the inverted red mask (gland)
	labels = measure.label(inverted_red, connectivity=1)
	
	# Step 2: Find components touching the boundary of the red mask (gland)
	rows, cols = red_mask.shape
	boundary_labels = np.unique(
		np.concatenate([
			labels[0, :], labels[rows - 1, :], labels[:, 0], labels[:, cols - 1]
		])
	)
	
	# Create a mask for components not touching the boundary of the red mask
	mask = np.isin(labels, boundary_labels, invert=True)
	
	# Apply the mask to the inverted red mask to get enclosed black regions
	enclosed_black = inverted_red & mask
	
	# Step 3: Ensure these regions do not touch any part of the blue mask
	# Dilate blue mask to ensure we catch touching regions
	blue_dilated = morphology.binary_dilation(blue_mask)
	
	# Remove black regions that touch the blue mask
	final_mask = enclosed_black & ~blue_dilated
	
	return final_mask

def edge_enclosed_black_regions(gland_mask, edge_mask):

	combined_mask = np.logical_or(gland_mask, edge_mask)

	# Step 1: Identify black regions in the (gland + edge)
	inverted_combined = ~combined_mask
	
	# Label connected components in the inverted red mask (gland)
	labels = measure.label(inverted_combined, connectivity=1)
	
	# Step 2: Find components touching the boundary of box (gland)
	rows, cols = combined_mask.shape
	boundary_labels = np.unique(
		np.concatenate([
			labels[0, :], labels[rows - 1, :], labels[:, 0], labels[:, cols - 1]
		])
	)
	
	# Create a mask for components not touching the boundary of the box
	enclosed_mask = np.isin(labels, boundary_labels, invert=True)
	
	# Apply the mask to the inverted red mask to get enclosed black regions
	enclosed_black = inverted_combined & enclosed_mask

	# Step 3: Find enclosed regions that touch the internal structure (gland)
	# Dilate the gland mask to make sure touching regions are detected
	gland_dilated = morphology.binary_dilation(gland_mask)

	# Find labels of regions touching the gland (internal structure)
	touching_labels = np.unique(labels[gland_dilated & enclosed_black])

	# Create a mask for all regions that touch the gland
	touching_internal_mask = np.isin(labels, touching_labels)

	return touching_internal_mask


def combine_components(main_image, additional_image):
	structure = ndi.generate_binary_structure(3, 1)
	combined_image = np.copy(main_image)
	
	# Find all components in the additional image that connect to the main image
	connected_mask = ndi.binary_dilation(main_image, structure=structure)
	combined_image |= additional_image & connected_mask
	
	return combined_image


def create_ellipsoidal_structuring_element(spacing, maxradius):
	# Calculate radius for each axis based on spacing, allowing zero where needed
	scaled_radius = [round(maxradius / s) for s in spacing]

	# Define the diameter for each dimension, setting a minimum of 1 for the shape where radius is zero
	diameter = [2 * r + 1 if r > 0 else 1 for r in scaled_radius]
	center = [d // 2 for d in diameter]

	# Create an empty grid and fill in the ellipsoidal structuring element
	ellipsoid = np.zeros(diameter, dtype=np.uint8)

	for x in range(diameter[0]):
		for y in range(diameter[1]):
			for z in range(diameter[2]):
				# Only calculate distance in dimensions with non-zero radius
				distance = sum(((i - c) / r) ** 2 for i, c, r in zip((x, y, z), center, scaled_radius) if r > 0)
				if distance <= 1:
					ellipsoid[x, y, z] = 1

	return ellipsoid

def single_step_binary_closing_3D(image, spacing, maxradius):

	# Calculate anisotropic structuring element based on spacing
	ellipsoid = create_ellipsoidal_structuring_element(spacing,maxradius)

   # Perform binary closing with the ellipsoidal structuring element
	closed_image = morphology.binary_closing(image, footprint=ellipsoid)

	return closed_image

def single_step_binary_opening_3D(image, spacing, maxradius):

	# Calculate anisotropic structuring element based on spacing
	ellipsoid = create_ellipsoidal_structuring_element(spacing,maxradius)

   # Perform binary closing with the ellipsoidal structuring element
	closed_image = morphology.binary_opening(image, footprint=ellipsoid)

	return closed_image

def single_step_grayscale_opening_3D(image, spacing, maxradius):

	# Calculate anisotropic structuring element based on spacing
	ellipsoid = create_ellipsoidal_structuring_element(spacing,maxradius)

   # Perform binary closing with the ellipsoidal structuring element
	closed_image = morphology.opening(image, footprint=ellipsoid)

	return closed_image

def single_step_grayscale_closing_3D(image, spacing, maxradius):

	# Calculate anisotropic structuring element based on spacing
	ellipsoid = create_ellipsoidal_structuring_element(spacing,maxradius)

   # Perform binary closing with the ellipsoidal structuring element
	closed_image = morphology.closing(image, footprint=ellipsoid)

	return closed_image

def single_step_grayscale_reconstruction_3D(image, spacing, maxradius):

	# Calculate anisotropic structuring element based on spacing
	ellipsoid = create_ellipsoidal_structuring_element(spacing,maxradius)

   # Perform binary closing with the ellipsoidal structuring element
	eroded_image = morphology.erosion(image, footprint=ellipsoid)

	reconstructed = morphology.reconstruction(eroded_image, image)

	return reconstructed

def single_step_binary_erosion_3D(image, spacing, maxradius):

	# Calculate anisotropic structuring element based on spacing
	ellipsoid = create_ellipsoidal_structuring_element(spacing,maxradius)

   # Perform binary closing with the ellipsoidal structuring element
	closed_image = morphology.binary_erosion(image, footprint=ellipsoid)

	return closed_image

def single_step_binary_dilation_3D(image, spacing, maxradius):

	# Calculate anisotropic structuring element based on spacing
	ellipsoid = create_ellipsoidal_structuring_element(spacing,maxradius)

   # Perform binary closing with the ellipsoidal structuring element
	closed_image = morphology.binary_dilation(image, footprint=ellipsoid)

	return closed_image

def remove_thin_connections(binary_mask, erosion_sizes, spacing, reference_mask=None, interior_points=None):
	"""
	Removes thin connections from a 3D binary mask through iterative erosion and
	watershed segmentation. Optionally uses a reference mask to retain only
	components overlapping with the reference mask.

	Parameters:
	- binary_mask: np.ndarray, 3D binary mask to process.
	- erosion_sizes: list of floats, sizes of erosions to apply iteratively.
	- spacing: tuple of floats, voxel spacing along each axis.
	- reference_mask: np.ndarray, optional, 3D binary mask from the previous frame
	  to retain overlapping components.

	Returns:
	- largest_region_mask: np.ndarray, binary mask of the largest connected region.
	"""      

	# Initialize the final labeled mask (starting empty)
	final_labels = np.zeros_like(binary_mask, dtype=int)

	if interior_points is not None:
		interior_mask = np.zeros_like(binary_mask, dtype=bool)
		if interior_points.size > 0:
			for z, y, x in interior_points.astype(int):
				if 0 <= z < binary_mask.shape[0] and 0 <= y < binary_mask.shape[1] and 0 <= x < binary_mask.shape[2]:
					interior_mask[z, y, x] = True

	# Step 1: Start with the largest erosion size
	for erosion_size in erosion_sizes:

		# Erode the original binary mask with the current erosion size
		eroded_mask = single_step_binary_erosion_3D(binary_mask, spacing, erosion_size)		
		label_counter = np.max(final_labels) + 1

		# Label the eroded regions
		new_labels = measure.label(eroded_mask,connectivity=3)

		# Step 2: Check each new label for overlaps with the existing final labels
		for label_id in range(1, np.max(new_labels) + 1):
			# Create a mask for the current new label
			current_label_mask = new_labels == label_id

			# If a reference_mask is provided, ensure overlapping regions are merged
			if reference_mask is not None:
				reference_overlap = np.any(reference_mask[current_label_mask])
				if reference_overlap:
					final_labels[current_label_mask] = 1  # Assign a unified label for all overlapping regions
					continue

			# If an interior_mask is provided, ensure overlapping regions are merged
			if interior_mask is not None:
				interior_overlap = np.any(interior_mask[current_label_mask])
				if interior_overlap:
					final_labels[current_label_mask] = 1  # Assign a unified label for all overlapping regions
					continue

			# Check if this label overlaps with any part of the existing final labels
			overlap = np.any(final_labels[current_label_mask] > 0)

			# If no overlap, add this new label to the final labels
			if not overlap:
				final_labels[current_label_mask] = label_counter
				label_counter +=1

	# Step 4: Perform watershed using the eroded labeled regions as seeds
	distance = ndi.distance_transform_edt(binary_mask)
	watershed_labels = segmentation.watershed(-distance, final_labels, mask=binary_mask)

	# Step 5: Extract the largest labeled area
	label_sizes = np.bincount(watershed_labels.ravel())
	largest_label = np.argmax(label_sizes[1:]) + 1  # Skip label 0 (background)
	largest_region_mask = watershed_labels == largest_label

	return largest_region_mask #watershed_labels

def fill_holes_in_mask(binary_mask):
	# Step 1: Invert the binary mask (holes become True in the inverted mask)
	inverted_mask = ~binary_mask

	# Step 2: Label connected components in the inverted mask
	labeled_inverted = measure.label(inverted_mask, connectivity=2)

	# Step 3: Generate a mask that is True at the boundary
	boundary_mask = np.zeros_like(binary_mask, dtype=bool)
	boundary_mask[0, :, :] = True   # Front face
	boundary_mask[-1, :, :] = True  # Back face
	boundary_mask[:, 0, :] = True   # Top face
	boundary_mask[:, -1, :] = True  # Bottom face
	boundary_mask[:, :, 0] = True   # Left face
	boundary_mask[:, :, -1] = True  # Right face

	# Step 4: Identify the labels that touch the boundary
	boundary_labels = np.unique(labeled_inverted[boundary_mask])

	# Step 5: Remove any labeled regions that overlap with the boundary
	hole_labels = np.setdiff1d(np.unique(labeled_inverted), boundary_labels)

	# Step 6: Create a mask of the holes and set them to True in the original mask
	hole_mask = np.isin(labeled_inverted, hole_labels)
	filled_mask = binary_mask | hole_mask  # Fill the holes by setting them to True

	return filled_mask


def watershed_2d_enclosed_space_ski_bg(raw_image, segmented_membrane, seed, bg, broaden=0):
	"""
	Perform watershed segmentation in 2D to find enclosed spaces within a membrane.
	
	Parameters:
	- raw_image: 2D numpy array, the original grayscale image.
	- segmented_membrane: 2D numpy array, a binary mask where the membrane is True (1) and background is False (0).
	- seed: region used as seed for watershed
	- bg: region that is background
	Returns:
	- labeled_image: 2D numpy array with labeled regions.
	"""

	# Invert the membrane to get the inside regions (membrane = 0, inside = 1).
	inverted_membrane = ~segmented_membrane

	markers = np.zeros(np.shape(raw_image),dtype=int)

	# Set the edges of the image as the outside markers (value 1)
	markers[bg] = 2
	markers[seed] = 1

	combined_image = np.copy(raw_image)
	combined_image[segmented_membrane] = 1

	# Perform watershed segmentation using the distance transform and markers.
	labels = segmentation.watershed(raw_image, markers=markers, mask=inverted_membrane, watershed_line=True)

	if broaden > 0:
		boundary = labels == 0
		boundary[segmented_membrane] = False
		broadened_boundary = morphology.binary_dilation(boundary, morphology.disk(broaden))  # Broaden by 1 pixels
		labels[broadened_boundary] = 0

	return labels


def rescale_image(image,bottom_cap,top_cap):
	Ifloor = np.percentile(image, 1.0)
	Icap = np.percentile(image, 99.5)

	image[image < bottom_cap*Icap] = bottom_cap*Icap
	image[image > top_cap*Icap] = top_cap*Icap
	image = (image - bottom_cap*Icap) 
	image = image / (Icap*(top_cap - bottom_cap))
	return image


def single_step_binary_closing_3D_itk(image_np, spacing, maxradius):
	"""
	Perform binary morphological closing on a 3D image, using physical radius scaling.
	
	Parameters:
		image_np (numpy.ndarray): Input binary 3D image (0 and 1).
		spacing (tuple): Tuple specifying voxel spacing (z, y, x).
		maxradius (float): Maximum physical radius for the structuring element.


	Returns:
		numpy.ndarray: Output binary 3D image after closing.
	"""
	Dimension = 3
	InputPixelType = itk.UC  # Unsigned Char
	ImageType = itk.Image[InputPixelType, Dimension]

	# Convert NumPy array to ITK image
	itk_image = itk.GetImageFromArray(image_np.astype(np.uint8))
	itk_image.SetSpacing(spacing)  # Set voxel spacing for anisotropic images

	# Create the ellipsoidal structuring element
	structuring_element = itk.FlatStructuringElement[Dimension].Ball(maxradius)

	# Set up the binary morphological closing filter
	closing_filter = itk.BinaryMorphologicalClosingImageFilter[ImageType, ImageType, itk.FlatStructuringElement[3]].New()
	closing_filter.SetInput(itk_image)
	closing_filter.SetKernel(structuring_element)
	closing_filter.SetForegroundValue(1) 

	# Execute the filter
	closing_filter.Update()
	closed_image_itk = closing_filter.GetOutput()

	# Convert back to NumPy
	closed_image_np = itk.GetArrayFromImage(closed_image_itk).astype(bool)
	return closed_image_np

def single_step_binary_opening_3D_itk(image_np, spacing, maxradius):
	"""
	Perform binary morphological opening on a 3D image, using physical radius scaling.
	
	Parameters:
		image_np (numpy.ndarray): Input binary 3D image (0 and 1).
		spacing (tuple): Tuple specifying voxel spacing (z, y, x).
		maxradius (float): Maximum physical radius for the structuring element.

	Returns:
		numpy.ndarray: Output binary 3D image after opening.
	"""
	Dimension = 3
	InputPixelType = itk.UC  # Unsigned Char (binary images are typically 0/255 in ITK)
	ImageType = itk.Image[InputPixelType, Dimension]

	# Convert NumPy array to ITK image
	itk_image = itk.GetImageFromArray(image_np.astype(np.uint8))
	itk_image.SetSpacing(spacing)  # Set voxel spacing for anisotropic images

	# Create the ellipsoidal structuring element
	structuring_element = itk.FlatStructuringElement[Dimension].Ball(maxradius)

	# Set up the binary morphological opening filter
	opening_filter = itk.BinaryMorphologicalOpeningImageFilter[ImageType, ImageType, itk.FlatStructuringElement[3]].New()
	opening_filter.SetInput(itk_image)
	opening_filter.SetKernel(structuring_element)
	opening_filter.SetForegroundValue(1)  # Specify the binary foreground value (assumes input is 0/1)

	# Execute the filter
	opening_filter.Update()
	opened_image_itk = opening_filter.GetOutput()

	# Convert back to NumPy
	opened_image_np = itk.GetArrayFromImage(opened_image_itk).astype(bool)
	return opened_image_np

def single_step_binary_dilation_3D_itk(image_np, spacing, maxradius):
	"""
	Perform binary morphological dilation on a 3D image, using physical radius scaling.
	
	Parameters:
		image_np (numpy.ndarray): Input binary 3D image (0 and 1).
		spacing (tuple): Tuple specifying voxel spacing (z, y, x).
		maxradius (float): Maximum physical radius for the structuring element.

	Returns:
		numpy.ndarray: Output binary 3D image after dilation.
	"""
	Dimension = 3
	InputPixelType = itk.UC  # Unsigned Char
	ImageType = itk.Image[InputPixelType, Dimension]

	# Convert NumPy array to ITK image
	itk_image = itk.GetImageFromArray(image_np.astype(np.uint8))
	itk_image.SetSpacing(spacing)  # Set voxel spacing for anisotropic images

	# Create the ellipsoidal structuring element
	structuring_element = itk.FlatStructuringElement[Dimension].Ball(maxradius)

	# Set up the binary morphological dilation filter
	dilation_filter = itk.BinaryDilateImageFilter[ImageType, ImageType, itk.FlatStructuringElement[3]].New()
	dilation_filter.SetInput(itk_image)
	dilation_filter.SetKernel(structuring_element)
	dilation_filter.SetForegroundValue(1)  # Specify the binary foreground value (assumes input is 0/1)

	# Execute the filter
	dilation_filter.Update()
	dilated_image_itk = dilation_filter.GetOutput()

	# Convert back to NumPy
	dilated_image_np = itk.GetArrayFromImage(dilated_image_itk).astype(bool)
	return dilated_image_np

def single_step_binary_erosion_3D_itk(image_np, spacing, maxradius):
	"""
	Perform binary morphological erosion on a 3D image, using physical radius scaling.
	
	Parameters:
		image_np (numpy.ndarray): Input binary 3D image (0 and 1).
		spacing (tuple): Tuple specifying voxel spacing (z, y, x).
		maxradius (float): Maximum physical radius for the structuring element.

	Returns:
		numpy.ndarray: Output binary 3D image after erosion.
	"""
	Dimension = 3
	InputPixelType = itk.UC  # Unsigned Char (binary images are typically 0/255 in ITK)
	ImageType = itk.Image[InputPixelType, Dimension]

	# Convert NumPy array to ITK image
	itk_image = itk.GetImageFromArray(image_np.astype(np.uint8))
	itk_image.SetSpacing(spacing)  # Set voxel spacing for anisotropic images

	# Create the ellipsoidal structuring element
	structuring_element = itk.FlatStructuringElement[Dimension].Ball(maxradius)

	# Set up the binary morphological erosion filter
	erosion_filter = itk.BinaryErodeImageFilter[ImageType, ImageType, itk.FlatStructuringElement[3]].New()
	erosion_filter.SetInput(itk_image)
	erosion_filter.SetKernel(structuring_element)
	erosion_filter.SetForegroundValue(1)  # Specify the binary foreground value (assumes input is 0/1)

	# Execute the filter
	erosion_filter.Update()
	eroded_image_itk = erosion_filter.GetOutput()

	# Convert back to NumPy
	eroded_image_np = itk.GetArrayFromImage(eroded_image_itk).astype(bool)
	return eroded_image_np

def PNA_based_preprocessing(image, image2, normalized_spacing, preprocess_image2=True):

	if preprocess_image2:
		# Find region of shared high intensity 
		t1 = threshold_otsu(image)
		t2 = threshold_otsu(image2)

		bin1 = image > t1
		bin2 = image2 > t2

		bin_c = np.logical_and(bin1,bin2)
		bin_o = largest_component(bin_c)
		bin_o =	single_step_binary_dilation_3D(bin_o, normalized_spacing, 2)
		distance = ndi.distance_transform_edt(~bin_o)
		x = 5
		feathered_mask = np.clip( (x - distance) / x, 0, 1)
		feathered_mask[bin_o] = 1  # Preserve the original mask

		# Preprocess image2
		image2_preprocess = np.clip(image2 - image, 0, None)
		image2_preprocess = image*feathered_mask + image2_preprocess*(1-feathered_mask)

	# Extract vasculature
		image_vasculature = np.clip(image - image2_preprocess, 0, None)
	else:
		image_vasculature = np.clip(image - image2, 0, None)

	# Extract the membrane
	image_membrane = np.clip(image - image_vasculature, 0, None)

	return image_membrane


def mask_based_preprocessing(image, image_PNA, mask, normalized_spacing, r_dilate=5, r_feather=10):

	image_vasculature = np.clip(image - image2, 0, None)
	image_membrane = np.clip(image - image_vasculature, 0, None)

	dilated_mask = single_step_binary_dilation_3D_itk(mask, normalized_spacing, r_dilate)

	distance = ndi.distance_transform_edt(~dilated_mask)
	feathered_mask = np.clip( (r_feather - distance) / r_feather, 0, 1)
	feathered_mask[dilated_mask] = 1  # Preserve the original mask			

	image_PNA = image_PNA*feathered_mask
	image_membrane = image_membrane*feathered_mask

	background_mask = np.invert(dilated_mask)

	return image_membrane, image_PNA, background_mask

def membrane_segmentation(plateness, global_threshold, axis, membrane_acceptance_level=100):
	"""
	Perform binarization using vesselness_downhill_binarization on slices along the specified axis.

	Parameters:
		data (numpy.ndarray): Output array to store the results.
		plateness (numpy.ndarray): Input plateness array.
		image_PNA (numpy.ndarray): Input image array.
		global_threshold (float): Threshold value.
		axis (int): Axis along which to process planes (0 for x, 1 for y, 2 for z).
	"""
	def process_single_plane(index):
		if axis == 0:
			return vesselness_downhill_binarization(plateness[index, :, :], global_threshold=global_threshold, membrane_acceptance_level=membrane_acceptance_level)
		elif axis == 1:
			return vesselness_downhill_binarization(plateness[:, index, :], global_threshold=global_threshold, membrane_acceptance_level=membrane_acceptance_level)
		elif axis == 2:
			return vesselness_downhill_binarization(plateness[:, :, index], global_threshold=global_threshold, membrane_acceptance_level=membrane_acceptance_level)

	data = np.zeros_like(plateness, dtype=bool)

	# Determine the dimension to iterate over
	dim = plateness.shape[axis]

	# Parallel processing of planes
	results = Parallel(n_jobs=-1)(
		delayed(process_single_plane)(index) for index in range(dim)
	)

	# Store the results in the appropriate axis
	for index, result in enumerate(results):
		if axis == 0:
			data[index, :, :] = result
		elif axis == 1:
			data[:, index, :] = result
		elif axis == 2:
			data[:, :, index] = result

	return data

def parallel_watershed_processing(image, membrane_mask, gland_mask, background_mask, axis):
	"""
	Perform parallel watershed segmentation for a specific axis.

	Parameters:
		image, membrane_mask, gland_mask, background_mask (numpy.ndarray): Input arrays.
		axis (int): Axis along which to process planes (0 for x-z, 1 for y-z).

	Returns:
		numpy.ndarray: Output data array after processing.
	"""
	def process_plane(index):
		if axis == 1:  # x-z slices
			slice_raw = image[:, index, :]
			slice_membrane = membrane_mask[:, index, :]
			slice_gland = gland_mask[:, index, :]
			slice_bg = background_mask[:, index, :]
		elif axis == 2:  # y-z slices
			slice_raw = image[:, :, index]
			slice_membrane = membrane_mask[:, :, index]
			slice_gland = gland_mask[:, :, index]
			slice_bg = background_mask[:, :, index]

		labels = watershed_2d_enclosed_space_ski_bg(slice_raw, slice_membrane, slice_gland, slice_bg)
		return index, labels

	# Determine output shape and dimension to iterate
	data_shape = image.shape
	dim = data_shape[axis]

	# Parallel processing of planes
	results = Parallel(n_jobs=-1)(
		delayed(process_plane)(index) for index in range(dim)
	)

	# Collect results into an output array
	data = np.zeros(data_shape, dtype=int)
	for index, labels in results:
		if axis == 1:
			data[:, index, :] = labels
		elif axis == 2:
			data[:, :, index] = labels

	return data

def extract_free_floating_components(inverted_image, axis):
	"""
	Extract free-floating components in parallel along the specified axis.

	Parameters:
		inverted_image (numpy.ndarray): Input binary image.
		axis (int): Axis along which to process planes (0 for x-z, 1 for y-z).

	Returns:
		numpy.ndarray: Boolean array with extracted free-floating components.
	"""
	def process_plane(index):
		if axis == 1:  # x-z slices
			return index, freefloating_components(inverted_image[:, index, :])
		elif axis == 2:  # y-z slices
			return index, freefloating_components(inverted_image[:, :, index])

	# Determine output shape and dimension to iterate
	output = np.zeros_like(inverted_image, dtype=bool)
	dim = inverted_image.shape[axis]

	# Parallel processing of planes
	results = Parallel(n_jobs=-1)(
		delayed(process_plane)(index) for index in range(dim)
	)

	# Collect results into the output array
	for index, component in results:
		if axis == 1:
			output[:, index, :] = component
		elif axis == 2:
			output[:, :, index] = component

	return output

def parallel_enclosed_regions_z(free_floating_complete, membrane_mask_coarse):
	"""
	Extract enclosed black regions in parallel along the z-axis.

	Parameters:
		free_floating_complete (numpy.ndarray): Input binary image.
		membrane_mask_coarse (numpy.ndarray): Mask image for processing.

	Returns:
		numpy.ndarray: Boolean array with processed slices.
	"""
	def process_plane(index):
		return index, enclosed_black_regions(free_floating_complete[index, :, :], membrane_mask_coarse[index, :, :])

	# Determine output shape
	output = np.zeros_like(free_floating_complete, dtype=bool)
	dim = free_floating_complete.shape[0]

	# Parallel processing of planes
	results = Parallel(n_jobs=-1)(
		delayed(process_plane)(index) for index in range(dim)
	)

	# Collect results into the output array
	for index, component in results:
		output[index, :, :] = component

	return output

def parallel_edge_enclosed_regions_z(free_floating_components, membrane_mask_coarse):
	"""
	Process edge-enclosed black regions in parallel along the z-axis.

	Parameters:
		free_floating_components (numpy.ndarray): Input binary image.
		membrane_mask_coarse (numpy.ndarray): Mask image for processing.

	Returns:
		numpy.ndarray: Boolean array with processed slices.
	"""
	def process_plane(index):
		return index, edge_enclosed_black_regions(free_floating_components[index, :, :], membrane_mask_coarse[index, :, :])

	# Determine output shape
	output = np.zeros_like(free_floating_components, dtype=bool)
	dim = free_floating_components.shape[0]

	# Parallel processing of planes
	results = Parallel(n_jobs=-1)(
		delayed(process_plane)(index) for index in range(dim)
	)

	# Collect results into the output array
	for index, component in results:
		output[index, :, :] = component

	return output

def parallel_watershed_processing_z(gradient, membrane_mask_coarse, gland_mask, background_mask, start_plane, end_plane):
	"""
	Perform watershed segmentation in parallel along the z-axis for specified planes.

	Parameters:
		gradient, membrane_mask_coarse, gland_mask, background_mask (numpy.ndarray): Input arrays.
		start_plane (int): Starting plane index.
		end_plane (int): Ending plane index.

	Returns:
		numpy.ndarray: Output array with watershed segmentation results.
	"""
	def process_plane(zplane):
		slice_raw = gradient[zplane, :, :]
		slice_membrane = membrane_mask_coarse[zplane, :, :]
		slice_gland = gland_mask[zplane, :, :]
		slice_bg = background_mask[zplane, :, :]
		labels = watershed_2d_enclosed_space_ski_bg(slice_raw, slice_membrane, slice_gland, slice_bg, broaden=1)
		return zplane, labels

	# Initialize output array
	output = np.zeros_like(gradient, dtype=int)
	output[gland_mask] = 1

	# Parallel processing of planes
	results = Parallel(n_jobs=-1)(
		delayed(process_plane)(zplane) for zplane in range(start_plane, end_plane)
	)

	# Collect results into the output array
	for zplane, labels in results:
		output[zplane, :, :] = labels

	return output

##################################
# Main loop
##################################

if __name__ == "__main__":

	FileFolder = "241104-SMG2"
	SaveFolder = "241104-SMG2"
	FileName = "2_SMG2_t"
	ImageName = FileName
	FilePath1 = f"./Data/{FileFolder}/{FileName}split_c1_rcan/"
	FilePath2 = f"./Data/{FileFolder}/{FileName}split_c2_rcan/"

	ResultPath = f"./Results/{SaveFolder}/"
	ResultDir = Path(ResultPath)
	ResultDir.mkdir(exist_ok=True)

	AnnotationPath = f"./Annotation/{SaveFolder}/"

	# Dimensions of the voxels
	voxel_size_x = 1.0 
	voxel_size_y = 0.8545552
	voxel_size_z = 0.8545552

	# Spacing used for plotting in Napari
	original_spacing = np.array([voxel_size_x, voxel_size_y, voxel_size_z])
	normalized_spacing = original_spacing / 1.0

	# Range for analysis			
	start = 0
	stop = 272
	count = 0
	step = 1
	prev_flag = False
	restart_flag = False

	zdim = 166
	xdim = 512
	ydim = 512

	# segmentation parameters
	membrane_acceptance_level = 10
	r_dilate = 5
	r_feather = 10
	r_erode = 5
	prev_thresh = 0.05
	preprocess_image2 = False 
	preprocess_background = True

	Ifloor = 0.00
	Icap = 1.00
	bottom_cap = 0.00
	top_cap = 1.00

	reference_path = f"{AnnotationPath}/reference_mask.tif"
	reference_mask = tiff.imread(reference_path).astype(bool)

	for i in range(start,stop+1,step):

		##########################################################
		# LOAD-FRAME
		##########################################################

		print(f"\nLoading frame {i}/{stop} ({count})")
		global_start_time = time.time()

		# Pre-allocate the image
		try:
			tif_data1 = tiff.imread(f"{FilePath1}{ImageName}{i:03d}.tif")					
			image = tif_data1[:,:,:] / (2**16 - 1) 

			tif_data2 = tiff.imread(f"{FilePath2}{ImageName}{i:03d}.tif")					
			image2 = tif_data2[:,:,:] / (2**16 - 1)

		except FileNotFoundError:
			print(f"Could not found files for frame {i:03d}")
			continue

		# dimensions of single images
		zdim = image.shape[0]
		xdim = image.shape[1]
		ydim = image.shape[2]

#		showNapari2Im(image,image2, normalized_spacing)

		##########################################################
		# LOAD ANNOTATION AND MASKS
		##########################################################

		if restart_flag:
			frame = i-1
			FilePath = f"./Results/{SaveFolder}/thresholded_epithelial_{frame:03d}.npz"
			prev_gland = np.load(FilePath)['data']
			prev_flag = True
			print(f"Restart using previous gland: {frame}")
			restart_flag = False

		rough_path = f"{AnnotationPath}/rough_mask_frame_000.tif"
		rough_mask = tiff.imread(rough_path).astype(bool)

		# @@@ Not very robust in case of restart. Requires restart from frame with interior points. NOT USED.
		interior_path = Path(f"{AnnotationPath}/interior_points_frame_{i:03d}.csv")
		if interior_path.exists():
			interior_points = np.loadtxt(f"{AnnotationPath}interior_points_frame_{i:03d}.csv", delimiter=",", skiprows=1)

		##########################################################
		# PRE-PROCESSING
		##########################################################

		print("Rescaling intensity")

		# Capping and rescaling of data
		image = rescale_image(image,bottom_cap,top_cap)
		image2 = rescale_image(image2,bottom_cap,top_cap)

		# Minor Gaussian blur
		image = filters.gaussian(image,sigma=1/normalized_spacing)
		image2 = filters.gaussian(image2,sigma=1/normalized_spacing)

		image_PNA = np.copy(image2)
		if prev_flag:
			image_membrane, image_PNA, background_mask = mask_based_preprocessing(image, image_PNA, prev_gland, normalized_spacing, r_dilate=r_dilate, r_feather=r_feather)
		else:
			image_membrane = PNA_based_preprocessing(image, image_PNA, normalized_spacing, preprocess_image2=preprocess_image2)
			background_mask = np.invert(rough_mask)

		##########################################################
		# FILTERING
		##########################################################

		print("Applying Hessian-based objectness filter")

		# Enhance platelike structures
		plateness = plateness_computation_spacing(image_membrane, normalized_spacing)
		plateness = plateness/2**16

		if not prev_flag:
			plateness[~rough_mask] = 0

		# Enhance high-contrast regions
		local_mean1 = filters.gaussian(plateness,sigma = 3/normalized_spacing, preserve_range=True)
		plateness = np.clip(plateness - local_mean1,0,None)

		##########################################################
		# MEMBRANE-SEGMENTATION
		##########################################################
		# Also includes any other features that light up through the marker

		print("Applying Membrane segmentation")

		# First binarization step through morphological reconstruction of 2D slices
		global_threshold = threshold_otsu(plateness)

		# Process y-z planes
		data_x = membrane_segmentation(plateness, global_threshold, axis=1, membrane_acceptance_level=membrane_acceptance_level)

		# Process x-z planes
		data_y = membrane_segmentation(plateness, global_threshold, axis=2, membrane_acceptance_level=membrane_acceptance_level)

		# Find parts of the basement membrane that appear in both data_x and data_y
		membrane_mask = np.logical_or(data_x,data_y)
		
		#showNapari(membrane_mask,normalized_spacing)
		membrane_mask = largest_component(membrane_mask)

		#showNapari2Im(membrane_mask,plateness,normalized_spacing)

		membrane_mask[reference_mask] = True

		# Create coarse segmentation through binary closing of the membrane image		
		membrane_mask_coarse = single_step_binary_closing_3D_itk(membrane_mask, normalized_spacing, 8)

		if prev_flag:
			reference_gland = np.copy(prev_gland)
			reference_gland = single_step_binary_erosion_3D(reference_gland, normalized_spacing, r_erode)
			reference_gland[membrane_mask] = False
			membrane_mask_coarse[reference_gland] = False

		##########################################################
		# SEGMENTATION OF LARGEST ENCLOSED STRUCTURE
		##########################################################
		# Should represent the epithelial cells
		# We segment based on connected components in the binary image
		# We do the segmentation slice-by-slice as this is more robust to holes in the membrane.

		print("Applying Gland extraction")
		# Take negative of the membrane image 
		inverted_image = ~membrane_mask_coarse

		# Extract the free floating components on a per slide basis and merge
		free_floating_components_x = extract_free_floating_components(inverted_image, axis=1)
		free_floating_components_y = extract_free_floating_components(inverted_image, axis=2)
		free_floating_components = np.logical_or(free_floating_components_x, free_floating_components_y)

		if prev_flag:
			# Only keep overlapping pieces with reference gland
			free_floating_components = overlapping_components(free_floating_components, reference_gland)
			free_floating_components[reference_gland] = True
		else:
			free_floating_components = overlapping_components_with_points(free_floating_components, interior_points)

		# Do a second round based on the merger of the freefloating parts and the membrane
		inner_plus_membrane = ~np.logical_or(free_floating_components, membrane_mask_coarse)
		free_floating_components_x = extract_free_floating_components(inner_plus_membrane, axis=1)
		free_floating_components_y = extract_free_floating_components(inner_plus_membrane, axis=2)
		loose_parts = np.logical_or(free_floating_components_x, free_floating_components_y)
		free_floating_complete = np.logical_or(loose_parts, free_floating_components)

		if prev_flag:
			free_floating_complete = overlapping_components(free_floating_complete, reference_gland)
		else:
			free_floating_complete = overlapping_components_with_points(free_floating_complete, interior_points)

		# Go through the z-plane to merge any holes in the gland
		free_floating_components_z = parallel_enclosed_regions_z(free_floating_complete, membrane_mask_coarse)
		free_floating_components = np.logical_or(free_floating_components_z, free_floating_complete)				

		if prev_flag:
			free_floating_components = overlapping_components(free_floating_components, reference_gland)
		else:
			free_floating_components = overlapping_components_with_points(free_floating_components, interior_points)

		# Go through the z-plane to merge any region enclosed by the gland + the membrane
		free_floating_components_z = parallel_edge_enclosed_regions_z(free_floating_components, membrane_mask_coarse)
		free_floating_components = np.logical_or(free_floating_components_z, free_floating_components)				

		# Prune structures with very thin connections
		if prev_flag:
			pruned = remove_thin_connections(free_floating_components, [4,3], normalized_spacing, reference_mask=reference_gland, interior_points=interior_points)
		else:
			pruned = remove_thin_connections(free_floating_components, [4,3], normalized_spacing, interior_points=interior_points)

		free_floating_components =  single_step_binary_opening_3D(pruned, normalized_spacing, 2)

		if prev_flag:
			gland_mask = overlapping_components(free_floating_components, reference_gland)
		else:
			gland_mask = overlapping_components_with_points(free_floating_components, interior_points)

		# Compute distance transform
		distance = ndi.distance_transform_edt(~membrane_mask_coarse)
		gradient = -distance + image + image2 # Consider using the processed data

		start_plane = 30
		end_plane = zdim - 30

		labels = parallel_watershed_processing_z(
			gradient,
			membrane_mask_coarse,
			gland_mask,
			background_mask,
			start_plane,
			end_plane,
		)

		gland_mask[labels == 1] = 1

		##########################################################
		# IMPROVE SALIVARY GLAND (watershed)
		##########################################################

		print("Improving gland segmentation")

		# Process x-z slices
		data_x = parallel_watershed_processing(
			image,
			membrane_mask,
			gland_mask,
			background_mask,
			1,
		)
		gland_x = data_x == 1

		# Process y-z slices
		data_y = parallel_watershed_processing(
			image,
			membrane_mask,
			gland_mask,
			background_mask,
			2,
		)
		gland_y = data_y == 1

		gland_combined = np.logical_and(gland_x,gland_y) # Why not or?
		glandwithmembrane = np.logical_or(gland_combined, membrane_mask)
		glandwithmembrane = fill_holes_in_mask(glandwithmembrane)
		gland_combined = np.copy(glandwithmembrane)
		gland_combined[membrane_mask] = False
		gland_combined = overlapping_components_with_points(gland_combined, interior_points)

		##########################################################
		# POST-PROCESSING (watershed)
		##########################################################
		# We keep processing limited here, as we can always refine the binary image later

		print("Post processing")
		# start_time = time.time()

		# Prune vasculature etc.
		glandwithmembrane = np.logical_or(gland_combined, membrane_mask)
		pruned = remove_thin_connections(glandwithmembrane, [4,3], normalized_spacing, interior_points=interior_points)
		gland_combined = np.copy(pruned)
		gland_combined[membrane_mask] = False

		gland_combined = fill_holes_in_mask(gland_combined)
		gland_combined =  single_step_binary_closing_3D(gland_combined, normalized_spacing, 1)
		membrane_mask[gland_combined] = False

		glandwithmembrane = np.logical_or(gland_combined, membrane_mask)
		glandwithmembrane =  single_step_binary_opening_3D(glandwithmembrane, normalized_spacing, 2)
		gland_combined = np.copy(glandwithmembrane)
		gland_combined[membrane_mask] = False
		gland_combined = largest_component(gland_combined)

		if not prev_flag:
			prev_flag = True
		else:
			dv = np.sum(np.logical_xor(gland_combined,prev_gland))
			frac = dv / (xdim*ydim*zdim)
			if frac > prev_thresh:
				prev_flag = False
			print(f"dv/V:\t{frac}\t{prev_flag}")

		prev_gland = np.copy(gland_combined)

		##########################################################
		# SAVE RESULTS
		##########################################################

		np.savez_compressed(f"./Results/{SaveFolder}/thresholded_membrane_{i:03d}", data=membrane_mask)
		np.savez_compressed(f"./Results/{SaveFolder}/thresholded_epithelial_{i:03d}", data=gland_combined)

		del membrane_mask
		del inner_plus_membrane
		del loose_parts
		del free_floating_complete
		del free_floating_components
		del free_floating_components_x
		del free_floating_components_y
		del free_floating_components_z

		global_end_time = time.time()
		elapsed_time = global_end_time - global_start_time
		print(f"frame time:\t{elapsed_time:.2f}")

		count +=1