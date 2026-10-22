"""
Script for viewing binary images.
Soft Math Lab (Harvard University)
Author: Justin Tauber

Collaboration with Shaohe Wang
"""

import numpy as np
import napari

def add_bounding_box(viewer, shape, scale):
	z, y, x = shape
	sz, sy, sx = scale

	corners = np.array([
		[0,   0,   0],
		[0,   0,   x* sx],
		[0,   y* sy, 0],
		[0,   y* sy, x* sx],
		[z* sz, 0,   0],
		[z* sz, 0,   x* sx],
		[z* sz, y* sy, 0],
		[z* sz, y* sy, x* sx],
	])

	edges = [
		[0, 1], [0, 2], [1, 3], [2, 3],  # bottom
		[4, 5], [4, 6], [5, 7], [6, 7],  # top
		[0, 4], [1, 5], [2, 6], [3, 7],  # verticals
	]

	lines = [corners[e] for e in edges]
	viewer.add_shapes(lines, shape_type='line', edge_color='yellow', edge_width=1)

def showNapari(image, spacing):
	aspect_xz = spacing[0]/spacing[1]

	# Create an empty viewer
	viewer = napari.Viewer(ndisplay=3)
	
	# Add the image with scaling
	layer = viewer.add_image(
		image,
		scale=(aspect_xz, 1, 1),
		rendering='attenuated_mip' 
	)
	
	viewer.layers.selection.active = layer

	# Add bounding box
	add_bounding_box(viewer, image.shape, (aspect_xz, 1, 1))


	napari.run()

if __name__ == "__main__":

	voxel_size_x = 1.0 
	voxel_size_y = 1.3811
	voxel_size_z =  1.3811
	original_spacing = np.array([voxel_size_x, voxel_size_y, voxel_size_z])

	Nframes = 283
			
	start = 165
	stop = start+1
	count = start
	step = 1

	for frame in range(start,stop,step):

		print("Loading File")
		FileName = f"./Results/250916-SMG2-BB94-washout_cp_model0907/thresholded_epithelial_{frame:03d}.npz"
		epithelialCells = np.load(FileName)['data']

		showNapari(epithelialCells,original_spacing)

		count +=1

