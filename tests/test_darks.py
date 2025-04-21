#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: test_darks.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import pathlib

import numpy


def test_median_dark(data_dir: pathlib.Path, tmp_path: pathlib.Path):
    """Test the creation of a median dark frame."""

    from ccd.bias import create_median_bias
    from ccd.darks import create_median_dark

    median_bias_file = tmp_path / "median_bias.fits"
    median_dark_file = tmp_path / "median_dark.fits"

    # Get the list of bias files
    bias_files = sorted(pathlib.Path(data_dir).glob("Bias*.fit*"))

    # Create the median bias frame
    median_bias = create_median_bias(bias_files, str(median_bias_file))
    assert median_bias_file.exists(), "The median bias file was not created."
    assert median_bias is not None, "The median bias frame is None."

    # Create the median dark.
    dark_files = sorted(pathlib.Path(data_dir).glob("Dark*.fit*"))
    median_dark = create_median_dark(
        dark_files,
        median_bias_file,
        str(median_dark_file),
    )

    assert median_dark_file.exists(), "The median dark file was not created."
    assert median_dark is not None, "The median dark frame is None."

    assert isinstance(median_dark, numpy.ndarray), "The median dark is not an array."
    assert median_dark.ndim == 2, "The median dark is not 2D."
    assert median_dark.shape == (4096, 4109), "The median dark has the wrong shape."
