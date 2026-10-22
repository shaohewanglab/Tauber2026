"""
Algorithm for computing area and volume from processed triangular surface meshes
Soft Math Lab (Harvard University)
Author: Justin Tauber

Collaboration with Shaohe Wang

"""
import numpy as np
import os
import trimesh

Folder = "cp_model0907_processed"
subFolder = "241104-SMG4_cp_model0907_v2"
folderPath = "/mnt/d/PD2_BuddingMorphogenesis/PaperData"

# Just looping over all possible frames in the folder
Nframes = 500
results = []

for i in range(Nframes):
	fileName = f"{folderPath}/Results_MeshVoxelization/{subFolder}/mesh_epithelial_{i:03d}.ply"	
	if os.path.exists(fileName):
		surfaces = trimesh.load(fileName, force='mesh')
		meshes = surfaces.split(only_watertight=False)
		if len(meshes) == 1:
			mesh = meshes[0]
		else:
			mesh = max(meshes, key=lambda m: m.volume if m.is_watertight else m.area)
		
		if mesh.is_watertight:
			volume = mesh.volume
		else:
			volume = np.nan 
		area = mesh.area

		time = i*5

		results.append((i, time, volume, area))

np.savetxt(f'{folderPath}/Results_MeshVoxelization/{subFolder}/volume_area_{subFolder}.txt', results, fmt='%d\t%d\t%.6f\t%.6f', header='frame\ttime_min\tvolume_um3\tsurface_area_um2', comments='')