"""
Code for converting binary images stored in compressed numpy arrays (.npz format) 
to tiff files.

Soft Math Lab (Harvard University)
Author: Justin Tauber

Collaboration with Shaohe Wang

"""
import numpy as np
import tifffile


if __name__ == "__main__":

	voxel_size_x = 1.0 
	voxel_size_y = 0.8545552
	voxel_size_z = 0.8545552

	Nframes = 218
			
	start = 0
	stop = Nframes
	count = start
	step = 1

	dataPath = "/mnt/d/PD2_BuddingMorphogenesis/PaperData/"
	folder = "241104-SMG2_cp_model0907"

	for frame in range(start,stop+1,step):

		print("Loading File")
		try:
			FileName = f"{dataPath}/Results/{folder}/thresholded_epithelial_{frame:03d}.npz"
			epithelialCells = np.load(FileName)['data']
		except FileNotFoundError:
			print(f"Could not find frame {frame}")
			continue

		OutName = f"{dataPath}/Results/{folder}/thresholded_epithelial_{frame:03d}.tif"

		int_array =  epithelialCells.astype("uint8") * 255
		tifffile.imwrite(OutName,
		                 int_array, 
		                 imagej=True,
		                 resolution=(1/voxel_size_y, 1/voxel_size_z),
						metadata={
						'spacing': voxel_size_x,
						'unit': 'um'
						})
							

		print(frame)

		count +=1