#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: ptc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ImageNormalize, ZScaleInterval

def calculate_gain(files):
    """This function:

    - Accept a list of files that you need to calculate the gain
      (two files should be enough, but what kind?).
    - Read the files and calculate the gain in e-/ADU.
    - Return the gain in e-/ADU.

    """
    #Loading flats
    flat1 = fits.getdata(files[0]).astype('f4')
    flat2 = fits.getdata(files[1]).astype('f4')
    # maybe need to trim ?
    #Finding mean
    mean1 = np.mean(flat1)
    mean2 = np.mean(flat2)
    #Finding variance
    variance = np.var(flat1-flat2)
    #finding gain
    gain = (mean1+mean2)/variance
    #retunring result
    return float(gain)


def calculate_readout_noise(files, gain):
    """This function:

    - Accept a list of files that you need to calculate the readout noise
      (two files should be enough, but what kind?).
    - Accept the gain in e-/ADU as gain. This should be the one you calculated
      in calculate_gain.
    - Read the files and calculate the readout noise in e-.
    - Return the readout noise in e-.

    """
    # Load the bias frames
    bias1 = fits.getdata(files[0]).astype('f4')
    bias2 = fits.getdata(files[1]).astype('f4')
    # Standard deviation
    std_diff = np.std(bias1 - bias2)
    # Readout noise
    readout_noise = (gain * std_diff) / np.sqrt(2)
    #return the result
    return readout_noise