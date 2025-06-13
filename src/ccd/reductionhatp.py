#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: reduction.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# i have hard coded into run reduction for the pytests. 

import glob
from astropy.io import fits
from astropy.table import vstack, Column
from pathlib import Path
from bias  import create_median_bias
from darks  import create_median_dark
from flats  import create_median_flat
from science import reduce_science_frame
from photometry import do_aperture_photometry, plot_radial_profile

def run_reduction(data_dir = None,
    out_dir = None,
    positions=[(511, 502), #my star 
               (458, 382), #comp star 1
              (389, 798) #comp star 2
              ],
    radii=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
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
#=========================================================================================================#
    # 1) Master bias
    bias_files = sorted(glob.glob(f"{data_dir}/Bias_BIN1_*.fits"))
    if not bias_files:
        print("⚠️  No bias frames found; skipping bias combination.")
    else:
        bias_master = f"{out_dir}/bias_master.fits"
        create_median_bias(bias_files, bias_master)
    print("bias done")
#=========================================================================================================#
    # 2) Master dark
    dark_files = sorted(glob.glob(f"{data_dir}/Dark_BIN1_*.fits"))
    if not dark_files:
        print("⚠️  No dark frames found. skipping")
    else:
        dark_master = f"{out_dir}/dark_master.fits"
        create_median_dark(dark_files, bias_master, dark_master)
    print("darks done")
#=========================================================================================================#
    # 3) Master flat (using z‐band flats)
    flat_files = sorted(glob.glob(f"{data_dir}/domeflat_z_*.fits"))
    if not flat_files:
        print("⚠️  No flats frames found; skipping.")
    else:
        flat_master = f"{out_dir}/flat_master.fits"
        create_median_flat(flat_files, bias_master, flat_master, dark_master)
    print("flats done")
#=========================================================================================================#
    # 4) Reduce science frames (will change the name of the stars)
    sci_files = sorted(glob.glob(f"{data_dir}/HAT-P-36 b_z_*.fits"))
    all_photometry = [] # save all into a table
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
        # 5) Photometry on science frames
        phot_table = do_aperture_photometry(
            image=reduced_name,
            positions=positions,
            radii=radii,
            sky_radius_in=sky_radius_in,
            sky_annulus_width=sky_annulus_width
        )

        # Extracting observation time from FITS header. I needed it for timer series analysis.
        with fits.open(sci) as hdul:
            header = hdul[0].header
            # I will try to get BJD-obs, fallback to less accurate. I did see it in all the frames i look at.
            time_obs = (
                header.get('BJD-OBS') or
                header.get('JD-OBS') or
                header.get('JD') or
                header.get('DATE-OBS')
            )

        # Adding filename and time as columns for later use. 
        phot_table['frame_filename'] = base
        phot_table['time_obs'] = time_obs

        all_photometry.append(phot_table)
        print(f"photometry + time for {base} saved")
      
        # 6) Plot radial profile
        plot_radial_profile(
            aperture_photometry_data=phot_table,
            output_filename=f"{out_dir}/{base}_radial_profile.png"
        )
        print("radial plot is done")

# Stacking all tables, sorted by time and saving (will change the name of the stars)
    if all_photometry:
        master_table = vstack(all_photometry)
        master_table.sort('time_obs')
        master_table.write(f"{out_dir}/HAT-P-36 b_photometry_master.ecsv", format="ascii.ecsv", overwrite=True)
        print(f"Combined photometry saved to {out_dir}/HAT-P-36 b_photometry_master.ecsv")
    else:
        print("No photometry tables created!")

        
    print("All done! Products are in", out_dir)
    return
