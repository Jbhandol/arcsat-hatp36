#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: darks.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np
def create_median_dark(dark_list, bias_filename, median_dark_filename):
    """This function:

    - Accept a list of dark file paths to combine as dark_list.
    - Accept a median bias frame filename as bias_filename (the one you created using
      create_median_bias).
    - Read all the images in dark_list and create a list of 2D numpy arrays.
    - Read the bias frame.
    - Subtract the bias frame from each dark image.
    - Divide each dark image by its exposure time so that you get the dark current
      per second. The exposure time can be found in the header of the FITS file.
    - Use a sigma clipping algorithm to combine all the bias-corrected dark frames
      using the median and removing outliers outside 3-sigma for each pixel.
    - Save the resulting dark frame to a FITS file with the name median_dark_filename.
    - Return the median dark frame as a 2D numpy array.

    """
    #Step Get data ()
    bias_data = fits.getdata(bias_filename).astype('f4')
    dark_frames = []
    #Step getting exposure time for each file and normalizing it. Also put them in float 32 arrays
    for filename in dark_list:
            dark = fits.open(filename)
            data = dark[0].data.astype('f4')
            exptime = dark[0].header['EXPTIME']
            dark_corrected = (data - bias_data) / exptime
            dark_frames.append(dark_corrected)
    dark_stack = np.array(dark_frames)
    #Step sigma clipping
    clipped = sigma_clip(dark_stack, sigma=3.0, axis=0, cenfunc='median')
    #Step average out
    median_dark = np.ma.mean(clipped, axis=0).data
    #Step Saving
    hdu = fits.PrimaryHDU(data=median_dark)
    hdu.writeto(median_dark_filename, overwrite=True)
    #Step returning the result when called out
    return median_dark

