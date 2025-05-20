#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: flats.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ImageNormalize, ZScaleInterval

def create_median_flat(
    flat_list,
    bias_filename,
    median_flat_filename,
    dark_filename,
):
    """This function:

    - Accept a list of flat file paths to combine as flat_list. Make sure all
      the flats are for the same filter.
    - Accept a median bias frame filename as bias_filename (the one you created using
      create_median_bias).
    - Read all the images in flat_list and create a list of 2D numpy arrays.
    - Read the bias frame.
    - Subtract the bias frame from each flat image.
    - Optionally you can pass a dark frame filename as dark_filename and subtract
      the dark frame from each flat image (remember to scale the dark frame by the
      exposure time of the flat frame).
    - Use a sigma clipping algorithm to combine all the bias-corrected flat frames
      using the median and removing outliers outside 3-sigma for each pixel.
    - Create a normalised flat divided by the median flat value.
    - Save the resulting median flat frame to a FITS file with the name
      median_flat_filename.
    - Return the normalised median flat frame as a 2D numpy array.

    """
    # Step initialize
    flat_frames = []
    bias_data = fits.getdata(bias_filename).astype('f4')
    if dark_filename is not None:
        dark_data = fits.getdata(dark_filename).astype('f4')
    else:
        dark_data = None
    # Step subtract bias and darks and append frames
    for file in flat_list:
        flat = fits.open(file)
        data = flat[0].data.astype('f4')
        exptime = float(flat[0].header['EXPTIME'])
        # Subtract bias
        bias_corrected = data - bias_data
        # Subtract scaled dark
        if dark_data is not None:
            darknbias_corrected = bias_corrected - dark_data * exptime
            flat_frames.append(darknbias_corrected)
        else:
            flat_frames.append(bias_corrected)
    # Step stack them
    flat_stack = np.array(flat_frames)
    # Step sigma clipping
    clipped = sigma_clip(flat_stack, sigma=3.0, axis=0, cenfunc='median')
    # Step combiing
    median_flat = np.ma.mean(clipped, axis=0).data
    # Step Nomralize
    median_flat /= np.median(median_flat)
    # Step Saving
    # Step Saving
    if median_flat_filename:  # Only save if filename is not None
        hdu = fits.PrimaryHDU(data=median_flat)
        hdu.writeto(median_flat_filename, overwrite=True)
    # Step retunring
    return median_flat


def plot_flat(
    median_flat_filename,
    ouput_filename="median_flat.png",
    profile_ouput_filename="median_flat_profile.png",
):
    """This function must:

    - Accept a normalised flat file path as median_flat_filename.
    - Read the flat file.
    - Plot the flat frame using matplotlib.imshow with reasonable vmin and vmax
      limits. Save the plot to the file specified by output_filename.
    - Take the median of the flat frame along the y-axis. You'll end up with a
      1D array.
    - Plot the 1D array using matplotlib.
    - Save the plot to the file specified by profile_output_filename.

    """
    # Load normalized flat
    flat = fits.getdata(median_flat_filename)

    # Plot the flat image
    norm = ImageNormalize(flat, interval=ZScaleInterval())
    plt.imshow(flat, cmap='gray', origin='lower', norm=norm)
    plt.colorbar()
    plt.title("Normalized Flat Field")
    # Save and close
    if ouput_filename:
        plt.savefig(ouput_filename)
    plt.close()

    #Plot the y-axis median profile
    profile = np.median(flat, axis=0)

    plt.plot(profile)
    plt.title("Median Flat Profile (Y-axis)")
    plt.xlabel("X pixel")
    plt.ylabel("Median Value")
    plt.grid(True)
    if profile_ouput_filename:
        plt.savefig(profile_ouput_filename)
    plt.close()
    return
