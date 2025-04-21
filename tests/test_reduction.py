#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: test_reduction.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations


def test_reduction():
    """Test the main reduction script. Just checks that it runs without errors."""

    from ccd.reduction import run_reduction

    # Run the reduction script
    run_reduction()

    assert True  # To make sure the test is not skipped.
