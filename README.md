# ASTR 480 Assignments Template

Template to create ASTR 480 assignments.

# ARCSAT HAT-P-36 b Photometry Pipeline

This repo contains all **code** needed to reproduce the calibration,
photometry, and light-curve analysis described in our ASTR 480 paper.
(The full raw FITS set is >1 GB and remains on UW JupyterHub; see path
below.)

## Quick-start

```bash
git clone https://github.com/Jbhandol/arcsat-hatp36.git
cd arcsat-hatp36
python run_reduction.py          # calls reduction.py internally

Outputs (reduction_output/)

bias_master.fits, dark_master.fits, flat_master.fits
master_table.ecsv – fluxes for target + two comparison stars
at every tested aperture
trimmed / normalised light-curve PNGs


For Further analysis

Open the Jupyter notebook
analysis_hat-p-36.ipynb
to load master_table.ecsv, pick the optimal aperture, and plot the
figures shown in the paper.

Full raw and calibration frames (≈800 MB) remain on UW JupyterHub at: https://jupyter.rttl.uw.edu/2025-spring-astr-480-a/hub/user-redirect/lab/tree/work/ccd-reductions-Jbhandol
