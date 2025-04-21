#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: test_ptc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import pathlib


def test_calculate_gain(data_dir: pathlib.Path):
    """Test the calculation of the gain in e-/ADU."""

    from ccd.ptc import calculate_gain

    files = [
        data_dir / "AutoFlat-PANoRot-r-Bin1-001.fit",
        data_dir / "AutoFlat-PANoRot-r-Bin1-002.fit",
    ]

    gain = calculate_gain(files)

    assert isinstance(gain, float), "The gain is not a float."
    assert gain > 0, "The gain is not positive."


def test_calculate_readout_noise(data_dir: pathlib.Path):
    """Test the calculation of the readout noise in e-."""

    from ccd.ptc import calculate_readout_noise

    files = [
        data_dir / "Bias-S001-R001-C001-NoFilt.fit",
        data_dir / "Bias-S001-R001-C002-NoFilt.fit",
    ]

    gain = 1.0  # Assume this is the gain value for the test purpose.

    readout_noise = calculate_readout_noise(files, gain)

    assert isinstance(readout_noise, float), "The readout noise is not a float."
    assert readout_noise > 0, "The readout noise is not positive."
