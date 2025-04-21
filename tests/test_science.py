#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: test_science.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import pathlib

import numpy


def test_reduce_science_frame(data_dir: pathlib.Path, tmp_path: pathlib.Path):
    """Test the creation of a median flat frame."""

    from ccd.bias import create_median_bias
    from ccd.darks import create_median_dark
    from ccd.flats import create_median_flat
    from ccd.science import reduce_science_frame

    median_bias_file = tmp_path / "median_bias.fits"
    median_dark_file = tmp_path / "median_dark.fits"
    median_flat_file = tmp_path / "median_flat.fits"
    science_reduced_file = tmp_path / "science_reduced.fits"

    # Get the lists of files
    bias_files = sorted(pathlib.Path(data_dir).glob("Bias*.fit*"))
    flat_files = sorted(pathlib.Path(data_dir).glob("AutoFlat*.fit*"))
    dark_files = sorted(pathlib.Path(data_dir).glob("Dark*.fit*"))

    # Create the median bias frame
    median_bias = create_median_bias(bias_files, str(median_bias_file))
    assert median_bias_file.exists(), "The median bias file was not created."
    assert median_bias is not None, "The median bias frame is None."

    # Create the median dark frame
    median_dark = create_median_dark(
        dark_files,
        median_bias_file,
        str(median_dark_file),
    )

    assert median_dark_file.exists(), "The median dark file was not created."
    assert median_dark is not None, "The median dark frame is None."

    # Create the median flat frame
    median_flat = create_median_flat(
        flat_files,
        str(median_bias_file),
        str(median_flat_file),
        dark_filename=str(median_dark_file),
    )

    assert median_flat_file.exists(), "The median flat file was not created."
    assert median_flat is not None, "The median flat frame is None."

    # Science frame to reduce
    science_files = pathlib.Path(data_dir) / "kelt-16-b-S001-R001-C084-r.fit"

    # Reduce the science frame
    data = reduce_science_frame(
        science_filename=str(science_files),
        median_bias_filename=str(median_bias_file),
        median_flat_filename=str(median_flat_file),
        median_dark_filename=str(median_dark_file),
        reduced_science_filename=str(science_reduced_file),
    )

    assert science_reduced_file.exists(), "The reduced science file was not created."

    assert data is not None, "The reduced science frame is None."
    assert isinstance(data, numpy.ndarray), "The reduced science frame is not an array."
    assert data.ndim == 2, "The reduced science frame is not a 2D array."
    assert data.shape == (4096, 4109), "The reduced science frame has the wrong shape."
