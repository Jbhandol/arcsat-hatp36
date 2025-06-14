from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable, hstack
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.aperture import ApertureStats
from photutils.centroids import centroid_com

def do_aperture_photometry(
    image,
    positions,
    radii,
    sky_radius_in,
    sky_annulus_width,
    centroiding=True,           # NEW: Turn on/off centroiding
    search_radius=10            # NEW: Radius (pixels) for centroid search box
):
    """
    this function is similar to what i had before. But i relaezide that there is shift in the flat fields. So i used the centroiding we discussed in the class. So this function performs aperture photometry on (optionally centroided) positions.

    We feed it:
        image: FITS filename (fully reduced science image).
        positions: List of (x, y) guess tuples for target and comparison stars.
        radii: List of aperture radii.
        sky_radius_in: Inner radius for sky annulus.
        sky_annulus_width: Width of sky annulus.
        centroiding: If True, refine positions using centroiding.
        search_radius: Pixel radius of search box for centroiding.

    Returns:
        Astropy table with background-subtracted flux for all stars, all apertures.
    """
    # Step 1: Load science image
    data = fits.getdata(image)
    refined_positions = []

    # Step 2: Centroid on each guess position (if enabled)
    for (x0, y0) in positions:
        if centroiding:
            # Extract small box
            cutout = data[
                int(y0)-search_radius:int(y0)+search_radius+1,
                int(x0)-search_radius:int(x0)+search_radius+1
            ]
            # Compute centroid
            cy, cx = centroid_com(cutout)
            x_star = x0 - search_radius + cx
            y_star = y0 - search_radius + cy
            refined_positions.append((x_star, y_star))
        else:
            refined_positions.append((x0, y0))

    # Step 3: Photometry (same as your original code)
    results = []
    for radius in radii:
        aperture = CircularAperture(refined_positions, r=radius)
        annulus = CircularAnnulus(refined_positions, r_in=sky_radius_in, r_out=sky_radius_in + sky_annulus_width)

        phot_ap = aperture_photometry(data, aperture)
        phot_an = aperture_photometry(data, annulus)

        annulus_stats = ApertureStats(data, annulus)
        bkg_median = annulus_stats.median
        bkg_total = bkg_median * aperture.area

        phot_ap[f'flux_r{radius}'] = phot_ap['aperture_sum'] - bkg_total

        results.append(phot_ap[[f'flux_r{radius}']])

    xs = [p[0] for p in refined_positions]
    ys = [p[1] for p in refined_positions]
    base_table = QTable([xs, ys], names=('xcenter', 'ycenter'))
    result = hstack([base_table] + results)
    result.meta['sky_radius_in'] = sky_radius_in
    return result

def plot_radial_profile(aperture_photometry_data, output_filename="radial_profile.png"):
    """This function must:

    - Accept a table of aperture photometry data as aperture_photometry_data. This
      is probably a photutils table, the result of do_aperture_photometry, but you
      can pass anything you like. The table/input data can contain a single target
      or multiple targets, but it must include multiple apertures.
    - Plot the radial profile of a target in the image using matplotlib. If you
      have multiple targets, label them correctly.
    - Plot a vertical line at the radius of the sky aperture used for the photometry.
    - Save the plot to the file specified in output_filename.

    """
    # Step 1: Identifing different radiuses from column names like 'flux_r5', 'flux_r10'
    radius_cols = [col for col in aperture_photometry_data.colnames if col.startswith('flux_r')]
    
    # Extract radius by splitting flux_r and 5. Mayeb not the most effeiceint way
    radii = [float(col.split('flux_r')[1]) for col in radius_cols]

    # Step 2: Plot one line per star (row)
    plt.figure()
    for i, row in enumerate(aperture_photometry_data):
        #get flux values for each radiuses/radii
        fluxes = [row[col] for col in radius_cols]
        plt.plot(radii, fluxes, marker='o', label=f"Star {i+1}")

    # grabing from meta  sky radius.
    sky_r = aperture_photometry_data.meta.get('sky_radius_in', None)
    # 4. If we have a sky radius, draw the line
    if sky_r is not None:
        plt.axvline(x=sky_r, color='gray', linestyle='--', label='Sky radius')

    plt.xlabel("Aperture Radius (pixels)")
    plt.ylabel("Background-subtracted Flux (ADU)")
    plt.title("Radial Profile of Aperture Photometry")
    plt.grid(True)
    plt.legend()
    plt.savefig(output_filename)
    plt.close()
    #pass