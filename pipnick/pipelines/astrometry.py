import logging
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs.wcs import FITSFixedWarning
from pathlib import Path
from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions import TimeoutError

from pipnick.utils.nickel_data import (bad_columns, ra_key, dec_key, fwhm_approx,
                                       ra_dec_units, fov_size_approx)
from pipnick.utils.dir_nav import organize_files

######################################
####  ASTROMETRY.NET CALIBRATION  ####
######################################

logger = logging.getLogger(__name__)

def astrometry_all(maindir, api_key, table_path=None, resolve=False,
                   excl_files=[], excl_objs=[], excl_filts=[]):
    """
    Runs astrometry.net calibration on all images input, and saves the calibrated
    fits files to astro_dir. Uses astrometry.net's API.

    Args:
        reddir (Path or str): path to directory containing images to calibrate
        resolve (bool): If True, re-solves images with previously generated local solves
    excl_files : list, optional
        List of file stems to exclude (exact match not necessary).
    excl_objs : list, optional
        List of object strings to exclude (exact match not necessary).
    excl_filts : list, optional
        List of filter names to exclude.
    
    Returns:
        list: list of relative paths (str) to all calibrated fits images
    """
    logger.info(f'---- astrometry_all() called on main directory {maindir}')
    ast = AstrometryNet()
    ast.api_key = api_key
    timeout = 60    # Seconds before timing out on initial solve
    
    maindir = Path(maindir)
    reddir = maindir / 'reduced'
    
    file_df = organize_files(reddir, table_path, 'astrometry',
                             excl_files, excl_objs, excl_filts)
    image_paths = file_df.paths
    
    # Modify images to remove the Nickel Telescope's bad columns
    mod_dir = maindir / 'processing' / 'mod-astro-input'  # Directory for modified images
    mod_dir.mkdir(parents=True, exist_ok=True)
    
    # Makes output folder if it doesn't already exist
    astro_dir = maindir / 'astrometric'
    wcs_dir = astro_dir / 'wcs'
    corr_dir = astro_dir / 'corr'
    astro_dir.mkdir(parents=True, exist_ok=True)
    wcs_dir.mkdir(parents=True, exist_ok=True)
    corr_dir.mkdir(parents=True, exist_ok=True)
    
    for obj_dir in reddir.iterdir():
        Path.mkdir(mod_dir / obj_dir.name, exist_ok=True)
        Path.mkdir(wcs_dir / obj_dir.name, exist_ok=True)
        Path.mkdir(corr_dir / obj_dir.name, exist_ok=True)
        
    mod_paths = make_cleaned_inputs(image_paths, mod_dir)

    def check_wcs_header(wcs_header):
        if wcs_header:
            # logger.info('')
            logger.info(f"Solution successful, saving WCS header to {output_path}")
            wcs_header.tofile(output_path, overwrite=True)
        else:
            # logger.warning('')
            logger.warning(f"Solution failed; skipping this image")

    # Ignore warnings about 'RADECSYS' header key being deprecated
    warnings.simplefilter("ignore", category=FITSFixedWarning)
    
    # Calibrate each image piece and collect output_paths
    calibrated_fits_paths = []
    unsolved_submission_ids = []
    for image_path in mod_paths:
        
        output_path_stem = image_path.stem.split('_')[0]
        obj_dir = image_path.parent.name
        output_path = wcs_dir / obj_dir / f"{output_path_stem}_astro.fits"
        if not resolve and output_path.exists():
            logger.info(f"Returning local copy of {image_path.name}'s solution; astrometry.net not used")
            calibrated_fits_paths.append(output_path)
            continue
        
        logger.info(f"Submitting image {image_path.name} to astrometry.net. This may take up to 60 seconds")
        try:
            wcs_header = ast.solve_from_image(image_path, fwhm=fwhm_approx,
                                              ra_key=ra_key, dec_key=dec_key,
                                              radius=fov_size_approx/60,  # Max error of RA/Dec guess, in deg
                                              scale_units='arcminwidth',
                                              scale_type='ul',
                                              scale_lower=fov_size_approx*0.8,
                                              scale_upper=fov_size_approx*1.2,
                                              solve_timeout=timeout,
                                              detect_threshold=8.0,
                                              publicly_visible='n',
                                              verbose=False)
            check_wcs_header(wcs_header)
        except TimeoutError as e:
            logger.info(f"No solution found on first try; will revisit after all images are submitted")
            unsolved_submission_ids.append(e.args[1])
        
    for submission_id in unsolved_submission_ids:
        try:
            wcs_header = ast.monitor_submission(submission_id,
                                                solve_timeout=120)
            check_wcs_header(wcs_header)
        except TimeoutError:
            logger.warning(f"No solution found on second try; skipping.")
    
    logger.info('---- astrometry_all() call ended')
    return calibrated_fits_paths


def make_cleaned_inputs(image_paths, mod_dir):
    # Modify images to remove the Nickel Telescope's bad columns
    logger.info("Zeroing out masked regions for faster astrometric solves")
    
    mod_paths = []
    for file in image_paths:
        mod_path = mod_dir / file.parent.name / (file.name.split('_')[0] + '_astroinput.fits')
        with fits.open(file) as hdul:
            data = np.array(hdul[0].data)
            try:
                # Creating modified FITS files w/ masked regions set to 0
                mask = np.array(hdul['MASK'].data)
                mask = mask.astype(bool)
                data[mask] = 0
            except KeyError:
                # If no mask in FITS file, sets bad columns = 0
                data[:, bad_columns] = 0
                logger.debug("No mask in FITS file--masking Nickel bad columns only")
            # Save the modified FITS file to the output directory
            fits.writeto(mod_path, data, hdul[0].header, overwrite=True)
        mod_paths.append(mod_path)
    
    return mod_paths


def get_astrometric_solves(image_paths, astro_dir, mode):
    """
    Returns any local copies of astrometric solves stored from previous runs of
    run_astrometry(). Skips if image has not yet been solved.

    Args:
        image_paths (list): list of relative paths to all images
        astro_dir (str): path to output image folder
        mode (str): Whether to return paths to calibrated image or .corr file w/ source table

    Returns:
        list: list of relative paths (str) to all calibrated fits images
    """
    
    logger.info("Returning local copies of astrometric solves; astrometry.net not used")
    calibrated_fits_paths = []
    astro_dir = Path(astro_dir)
    image_paths = [Path(image_path) for image_path in image_paths]
    for image_path in image_paths:
        output_path_stem = image_path.stem
        obj_dir = output_path.parent
        if mode == 'image':
            output_path = astro_dir / obj_dir.name / f"{output_path_stem[:-5]}.fits"
        elif mode == 'corr':
            output_path = astro_dir / obj_dir.name / f"{output_path_stem[:-5]}_corr.fits"
        if output_path.exists():
            logger.debug(f"Found calibrated image {Path(image_path).name}; appending to list")
            calibrated_fits_paths.append(output_path)
    return calibrated_fits_paths
