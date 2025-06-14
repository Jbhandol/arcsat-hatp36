#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: bias.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np
def create_median_bias(bias_list, median_bias_filename):
    """This function:

    - Accepts a list of bias file paths as bias_list.
    - Read each bias file and create a list of 2D numpy arrays.
    - Use a sigma clipping algorithm to combine all the bias frames using
      the median and removing outliers outside 3-sigma for each pixel.
    - Save the resulting median bias frame to a FITS file with the name
      median_bias_filename.
    - Return the median bias frame as a 2D numpy array.

    """
    #Step: Load data (here i will read each relevant FITS file, convert to float 32)
    bias_frames = []
    for filename in bias_list:
        data = fits.getdata(filename).astype('f4')  # convert to float32
        bias_frames.append(data)
    #Step: Make a 3D array (Numpy arrary)
    bias_3d_array = np.array(bias_frames)
    #Step: CLipping (use sigma_clip with median and sigma of 2 or 3)
    sigma_clipped_array = sigma_clip(bias_3d_array, sigma=2.5, axis=0, cenfunc='median')
    #Step: Averagingm (mean of unmasked data)
    median_bias = np.ma.mean(sigma_clipped_array, axis=0).data
    #Step: Saving (median_bias_filename)
    primary = fits.PrimaryHDU(data=median_bias, header=fits.Header())
    hdul = fits.HDUList([primary])
    hdul.writeto(median_bias_filename, overwrite=True)
    #Step: Return reuslt (shoudl return 2d array when called uponm)
    return median_bias