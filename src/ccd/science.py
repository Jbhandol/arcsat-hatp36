#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: science.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ImageNormalize, ZScaleInterval

def reduce_science_frame(
    science_filename,
    median_bias_filename,
    median_flat_filename,
    median_dark_filename,
    reduced_science_filename="reduced_science.fits",
):
    """This function must:

    - Accept a science frame filename as science_filename.
    - Accept a median bias frame filename as median_bias_filename (the one you created
      using create_median_bias).
    - Accept a median flat frame filename as median_flat_filename (the one you created
      using create_median_flat).
    - Accept a median dark frame filename as median_dark_filename (the one you created
      using create_median_dark).
    - Read all files.
    - Subtract the bias frame from the science frame.
    - Subtract the dark frame from the science frame. Remember to multiply the
      dark frame by the exposure time of the science frame. The exposure time can
      be found in the header of the FITS file.
    - Correct the science frame using the flat frame.
    - Optionally, remove cosmic rays.
    - Save the resulting reduced science frame to a FITS file with the filename
      reduced_science_filename.
    - Return the reduced science frame as a 2D numpy array.

    """
    # Step: Loading all the relevant calibration files
    bias = fits.getdata(median_bias_filename).astype('f4')
    flat = fits.getdata(median_flat_filename).astype('f4')
    dark = fits.getdata(median_dark_filename).astype('f4')
    # Step: Getting the science pic
    sci = fits.open(science_filename)
    sci_data = sci[0].data.astype('f4')
    exptime = float(sci[0].header['EXPTIME'])
    # Step Subtracting the bias
    reduced = sci_data - bias
    # Step Subtracting the dark
    reduced -= dark * exptime
    # Step Trimming to correct dimensions
    reduced = reduced
    #flat = flat[1000:3000, 1000:3000]
    # Step divide by normalized flat
    reduced /= flat
    # Step Save
    hdu = fits.PrimaryHDU(reduced)
    hdu.writeto(reduced_science_filename, overwrite=True)
    # Step remove cosmic rays maybe
    # Step return
    return reduced