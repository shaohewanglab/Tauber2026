"""
Convert segmented volumes saved as binary images into surface meshes
Used as input for voxelization protocol.

Soft Math Lab (Harvard University)
Author: Justin Tauber

Collaboration with Shaohe Wang
"""

import numpy as np
from scipy import linalg
from skimage import measure
from skimage import filters

import trimesh
import open3d as o3d

################################################
# Functions for mesh processing
################################################

def pad_with(vector, pad_width, iaxis, kwargs):
	pad_value = kwargs.get('padder', 0.0)
	vector[:pad_width[0]] = pad_value
	vector[-pad_width[1]:] = pad_value

def taubin_smoothing(mesh, iterations=1, lambda_=0.5, mu=-0.53):
	"""Apply Taubin smoothing to a mesh."""
	for it in range(iterations):
		vertices = mesh.vertices.copy()
		
		# First pass (Laplacian smoothing)
		laplacian = np.zeros_like(vertices)
		for i, v in enumerate(vertices):
			neighbors = mesh.vertex_neighbors[i]
			laplacian[i] = np.mean(vertices[neighbors], axis=0) - v
		mesh.vertices += lambda_ * laplacian
		
		# Second pass (inflation)
		vertices = mesh.vertices.copy()
		inflation = np.zeros_like(vertices)
		for i, v in enumerate(vertices):
			neighbors = mesh.vertex_neighbors[i]
			inflation[i] = np.mean(vertices[neighbors], axis=0) - v
		mesh.vertices += mu * inflation

	return mesh

# Convert trimesh mesh to open3d mesh
def trimesh_to_open3d(mesh):
	vertices = np.asarray(mesh.vertices)
	triangles = np.asarray(mesh.faces)
	o3d_mesh = o3d.geometry.TriangleMesh()
	o3d_mesh.vertices = o3d.utility.Vector3dVector(vertices)
	o3d_mesh.triangles = o3d.utility.Vector3iVector(triangles)
	return o3d_mesh

# Convert back from open3d mesh to trimesh
def open3d_to_trimesh(o3d_mesh):
	vertices = np.asarray(o3d_mesh.vertices)
	triangles = np.asarray(o3d_mesh.triangles)
	return trimesh.Trimesh(vertices=vertices, faces=triangles)


def skimage_to_open3d(verts, faces):
	mesh = o3d.geometry.TriangleMesh()   
	mesh.vertices = o3d.utility.Vector3dVector(verts)
	mesh.triangles = o3d.utility.Vector3iVector(faces)
	mesh.compute_vertex_normals()
	return mesh

if __name__ == "__main__":

	voxel_size_x = 1.0
	voxel_size_y = 1.3811
	voxel_size_z = 1.3811

	original_spacing = np.array([voxel_size_x, voxel_size_y, voxel_size_z])
	normalized_spacing = original_spacing / 1.0
	
	Nframes = 218
			
	start = 0
	stop = Nframes
	count = start
	step = 1
	rotate_mesh = False

	folder = "250916-SMG8-BB94-washout_cp_model0907"

	for frame in range(start,stop+1,step):

		print("Loading File")
		try:
			FileName = f"./Results/{folder}/thresholded_epithelial_{frame:03d}.npz"
			epithelialCells = np.load(FileName)['data']
		except FileNotFoundError:
			print(f"Could not find frame {frame}")
			continue
		bg = np.invert(epithelialCells)

		# Use a gaussian filter to create smoother edges 
		sigma = 0.5  # Adjust sigma to control the amount of smoothing
		epithelialCells = filters.gaussian(epithelialCells, sigma= 1.0 / normalized_spacing)
		bg = filters.gaussian(bg, sigma= 1.0 / normalized_spacing) 
		smoothed_image = np.clip(epithelialCells - bg, 0, None)

		# Create a surface mesh with marching cubes (slightly higher threshold to prevent merging of sharp features)
		padded = np.pad(epithelialCells, 1, pad_with)
		verts, faces, normals, values = measure.marching_cubes(padded,0.63,spacing = normalized_spacing, gradient_direction="ascent")
		verts = np.asarray(verts, dtype=np.float64)
		faces = np.asarray(faces, dtype=np.int32)
		  
		# Extract Largest Component
		o3d_mesh = skimage_to_open3d(verts, faces)
		labels = np.array(o3d_mesh.cluster_connected_triangles()[0])
		largest_component = np.argmax(np.bincount(labels))
		triangles_to_remove = [i for i, label in enumerate(labels) if label != largest_component]
		o3d_mesh.remove_triangles_by_index(triangles_to_remove)

		# Apply initial Taubin smoothing to remove high frequency noise
		smoothing_iterations = 10 
		taubin_lambda = 0.5 
		taubin_mu = -0.53 
		smoothed_mesh = o3d_mesh.filter_smooth_taubin(number_of_iterations=smoothing_iterations, lambda_filter=taubin_lambda, mu=taubin_mu)

		# Iteratively decimate the surface meshes followed by Taubin smoothing to reduce the number of points in the surface mesh
		n_iterations = 11
		decimate_factor = 0.75
		smoothing_iterations = 5

		for i in range(n_iterations):
			# Decimate Mesh with Quadric Decimation
			target_number_of_triangles = int(decimate_factor*np.asarray(smoothed_mesh.triangles).shape[0])  # Adjust as needed
			decimated_mesh = smoothed_mesh.simplify_quadric_decimation(target_number_of_triangles)          

			# Apply Taubin Smoothing
			smoothed_mesh = decimated_mesh.filter_smooth_taubin(number_of_iterations=smoothing_iterations, lambda_filter=taubin_lambda, mu=taubin_mu)
			print(f"smooth {i}")

		# Check if the mesh is watertight
		smoothed_mesh = open3d_to_trimesh(smoothed_mesh)
		print(smoothed_mesh.is_watertight)

		# Calculate genus
		Nvertices = len(smoothed_mesh.vertices)
		Nedges = len(smoothed_mesh.edges)/2
		Nfaces = len(smoothed_mesh.faces)
			
		chi = (Nvertices - Nedges + Nfaces)
		g = (2.-chi)/2

		print(g)

		# Save the mesh
		output_file = f'./Results/{folder}/mesh_epithelial_{frame:03d}.ply'
		smoothed_mesh.export(output_file)

		print(frame)

		count +=1

