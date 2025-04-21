#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: test_bias.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import pathlib

import numpy


def test_median_bias(data_dir: pathlib.Path, tmp_path: pathlib.Path):
    """Test the creation of a median bias frame."""

    from ccd.bias import create_median_bias

    test_file = tmp_path / "test_median_bias.fits"

    # Get the list of bias files
    bias_files = sorted(pathlib.Path(data_dir).glob("Bias*.fit*"))

    # Create the median bias frame
    median_bias = create_median_bias(bias_files, str(test_file))

    assert test_file.exists(), "The median bias file was not created."

    assert isinstance(median_bias, numpy.ndarray), "The median bias is not an array."
    assert median_bias.ndim == 2, "The median bias is not a 2D array."
    assert median_bias.shape == (4096, 4109), "The median bias has the wrong shape."
