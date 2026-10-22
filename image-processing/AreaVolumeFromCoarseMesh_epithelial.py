"""
Algorithm for computing area and volume from coarse triangular surface meshes
Soft Math Lab (Harvard University)
Author: Justin Tauber

Collaboration with Shaohe Wang

"""
import numpy as np
import os

import trimesh

# Parameters
Folder = "Results"
subFolders = ["250916-SMG1-BB94-washout","250916-SMG2-BB94-washout","250916-SMG3-BB94-washout","250916-SMG4-BB94-washout","250916-SMG5-BB94-washout","250916-SMG6-BB94-washout","250916-SMG7-BB94-washout","250916-SMG8-BB94-washout"]#,
folderPath = "/mnt/d/PD2_BuddingMorphogenesis/PaperData"

Nframes_array = [218,160,200,218,218,218,218,218]

for j, subFolder in enumerate(subFolders):

	subFolder = f"{subFolder}_cp_model0907"

	Nframes = Nframes_array[j]
	results = []

	for i in range(Nframes):
		fileName = f"{folderPath}/{Folder}/{subFolder}/mesh_epithelial_{i:03d}.ply"	
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

	np.savetxt(f'{folderPath}/{Folder}/{subFolder}/volume_area_{subFolder}.txt', results, fmt='%d\t%d\t%.6f\t%.6f', header='frame\ttime_min\tvolume_um3\tsurface_area_um2', comments='')
