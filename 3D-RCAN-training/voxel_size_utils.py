# MIT License
#
# Copyright (c) 2026 Shaohe Wang
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Read and write the voxel size of ImageJ TIFF stacks.

Used by the modified 3D-RCAN apply.py so that restored images keep the voxel
size and unit of their input. Missing calibration falls back to 1.0 (pixel
units) with a warning instead of failing.
"""

import numpy as np
import tifffile


def get_imagej_metadata(file_path):
    """
    Extract voxel size and unit information from an ImageJ format TIFF file.

    Parameters:
        file_path (str): Path to the TIFF file.

    Returns:
        dict: voxel_size_x, voxel_size_y, voxel_size_z (floats) and unit (str).
        Values missing from the file are set to 1.0 and 'pixel', with a warning.
    """
    metadata = {
        'voxel_size_x': None,
        'voxel_size_y': None,
        'voxel_size_z': None,
        'unit': None
    }

    with tifffile.TiffFile(file_path) as tif:
        # Check for ImageJ metadata
        imagej_metadata = tif.imagej_metadata
        if imagej_metadata:
            if 'spacing' in imagej_metadata:
                metadata['voxel_size_z'] = imagej_metadata['spacing']  # Z voxel size
            if 'unit' in imagej_metadata:
                metadata['unit'] = imagej_metadata['unit']  # Unit of measurement

        # Get XY resolutions from tags
        page = tif.pages[0]
        if page.tags.get('XResolution') and page.tags.get('YResolution'):
            # In TIFF files, XResolution (and YResolution) is stored as a
            # rational number (numerator, denominator) giving pixels per unit.
            x_res = page.tags['XResolution'].value
            y_res = page.tags['YResolution'].value
            # voxel sizes are in units per pixel
            if x_res[0] > 0 and y_res[0] > 0:
                metadata['voxel_size_x'] = 1.0 * x_res[1] / x_res[0]
                metadata['voxel_size_y'] = 1.0 * y_res[1] / y_res[0]

    missing = [k for k in ('voxel_size_x', 'voxel_size_y', 'voxel_size_z')
               if metadata[k] is None]
    if missing:
        print('Warning: no calibration for %s in %s; using 1.0' % (', '.join(missing), file_path))
        for k in missing:
            metadata[k] = 1.0
    if metadata['unit'] is None:
        metadata['unit'] = 'pixel'
    return metadata


def save_imagej_hyperstack(filename, image, metadata=None):
    """
    Save a (Z, Y, X) or (C, Z, Y, X) array as an ImageJ hyperstack with the
    voxel size and unit in metadata (as returned by get_imagej_metadata).
    Without metadata the voxel size is 1.0 pixel.
    """
    assert image.ndim in [3, 4]
    if metadata is None:
        metadata = {'voxel_size_x': 1.0, 'voxel_size_y': 1.0,
                    'voxel_size_z': 1.0, 'unit': 'pixel'}
    if image.ndim == 4:
        image = np.transpose(image, (1, 0, 2, 3))  # (C, Z, Y, X) -> (Z, C, Y, X)
    tifffile.imwrite(str(filename), image, imagej=True,
                     resolution=(1.0 / metadata['voxel_size_x'],
                                 1.0 / metadata['voxel_size_y']),
                     metadata={
                         'spacing': metadata['voxel_size_z'],
                         'unit': metadata['unit'],
                         'axes': 'ZYX' if image.ndim == 3 else 'ZCYX',
                     })
