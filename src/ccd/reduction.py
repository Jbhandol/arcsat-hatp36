#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: reduction.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# i have hard coded into run reduction for the pytests. 

import glob
from astropy.io import fits
from pathlib import Path
from ccd.bias  import create_median_bias
from ccd.darks  import create_median_dark
from ccd.flats  import create_median_flat
from ccd.science import reduce_science_frame
from ccd.photometry import do_aperture_photometry, plot_radial_profile

def run_reduction(data_dir = None,
    out_dir = None,
    positions=[(1305.57, 1022.5), (786.009, 1115.66)],
    radii=[5, 10, 15],
    sky_radius_in=20,
    sky_annulus_width=5):
    """This function must run the entire CCD reduction process. You can implement it
    in any way that you want but it must perform a valid reduction for the two
    science frames in the dataset using the functions that you have implemented in
    this module. Then perform aperture photometry on at least one of the science
    frames, using apertures and sky annuli that make sense for the data.

    No specific output is required but make sure the function prints/saves all the
    relevant information to the screen or to a file, and that any plots are saved to
    PNG or PDF files.

    """
    print("starting reduction")
    #finding correct dirs
    cwd = Path.cwd()
    if data_dir:
        # if someone passes a path in, use it
        data_dir = Path(data_dir)
    elif (cwd / "ccd_reductions_data").exists():
        # pytest and Classroom put data here
        data_dir = cwd / "ccd_reductions_data"
    else:
        # my local JupyterHub setup uses data/ccd_reductions_data. i should have put it in same folder. meh
        data_dir = cwd / "data" / "ccd_reductions_data"
    # outdir cleanup
    out_dir = Path(out_dir or cwd / "reduction_output")
    out_dir.mkdir(exist_ok=True)

    #doign actual reduction
    # 1) Master bias
    
    bias_files = sorted(glob.glob(f"{data_dir}/Bias-*.fit"))
    if not bias_files:
        print("⚠️  No bias frames found; skipping bias combination.")
    else:
        bias_master = f"{out_dir}/bias_master.fits"
        create_median_bias(bias_files, bias_master)
    print("bias done")
    # 2) Master dark
    dark_files = sorted(glob.glob(f"{data_dir}/Dark-*.fit"))
    if not dark_files:
        print("⚠️  No dark frames found. skipping")
    else:
        dark_master = f"{out_dir}/dark_master.fits"
        create_median_dark(dark_files, bias_master, dark_master)
    print("darks done")
    # 3) Master flat (using r‐band flats)
    flat_files = sorted(glob.glob(f"{data_dir}/AutoFlat-PANoRot-*-Bin1-*.fit"))
    if not flat_files:
        print("⚠️  No flats frames found; skipping.")
    else:
        flat_master = f"{out_dir}/flat_master.fits"
        create_median_flat(flat_files, bias_master, flat_master, dark_master)
    print("flats done")
    # 4) Reduce both science frames
    sci_files = sorted(glob.glob(f"{data_dir}/kelt-16-b-*.fit"))
    for sci in sci_files:
        base = sci.split("/")[-1].replace(".fit","")
        reduced_name = f"{out_dir}/{base}_reduced.fits"
        reduce_science_frame(
            science_filename=sci,
            median_bias_filename=bias_master,
            median_flat_filename=flat_master,
            median_dark_filename=dark_master,
            reduced_science_filename=reduced_name
        )
        print("sceince frame reduction done")
        # 5) Photometry on first science frame only (for now)
        phot_table = do_aperture_photometry(
            image=reduced_name,
            positions=positions,
            radii=radii,
            sky_radius_in=sky_radius_in,
            sky_annulus_width=sky_annulus_width
        )
        print("photometry done")
        # 6) Plot radial profile
        plot_radial_profile(
            aperture_photometry_data=phot_table,
            output_filename=f"{out_dir}/{base}_radial_profile.png"
        )
        print("radial plot is done")
    print("All done! Products are in", out_dir)
    return
