#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: flats.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


def create_median_flat(
    flat_list,
    bias_filename,
    median_flat_filename,
    dark_filename=None,
):
    """This function must:

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

    # This is a placeholder for the actual implementation.
    median_flat = None

    # See code in create_median_bias for how to create a new FITS file
    # from the resulting median flat frame.

    return median_flat


def plot_flat(median_flat_filename):
    """This function must:

    - Accept a normalised flat file path as median_flat_filename.
    - Read the flat file.
    - Plot the flat frame using matplotlib.imshow with reasonable vmin and vmax
      limits. Save the plot to a file called median_flat.png
    - Take the median of the flat frame along the y-axis. You'll end up with a
      1D array.
    - Plot the 1D array using matplotlib.
    - Save the plot to a file with the name median_flat_profile.png.

    """

    return
