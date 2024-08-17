import logging
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs.wcs import FITSFixedWarning
from pathlib import Path
from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions import TimeoutError

from pipnick.utils.nickel_data import (bad_columns, ra_key, dec_key, fwhm_approx,
                                       fov_size_approx)
from pipnick.utils.dir_nav import organize_files

######################################
####  ASTROMETRY.NET CALIBRATION  ####
######################################

logger = logging.getLogger(__name__)

def astrometry_all(maindir, api_key, table_path=None, resolve=False,
                   excl_files=[], excl_objs=[], excl_filts=[]):
    """
    Run astrometry.net calibration on all images in the specified directory.

    This function submits images to the astrometry.net API for calibration, and 
    saves the WCS header results in a dedicated directory. It handles the 
    Nickel Telescope's bad columns and optionally resolves previously solved images.

    Parameters
    ----------
    maindir : str or Path
        Path to the main directory containing images to calibrate.
    api_key : str
        API key for astrometry.net.
    table_path : str or Path, optional
        Path to a table containing metadata for the images.
    resolve : bool, optional
        If True, re-solves images with previously generated local solutions.
    excl_files : list of str, optional
        List of file stems to exclude (exact match not necessary).
    excl_objs : list of str, optional
        List of object strings to exclude (exact match not necessary).
    excl_filts : list of str, optional
        List of filter names to exclude.

    Returns
    -------
    list of Path
        List of relative paths to all calibrated FITS images.
    """
    logger.info(f'---- astrometry_all() called on main directory {maindir}')
    ast = AstrometryNet()
    ast.api_key = api_key
    timeout = 60  # Seconds before timing out on initial solve
    
    maindir = Path(maindir)
    reddir = maindir / 'reduced'
    
    # Organize files based on exclusion criteria and directory structure
    file_df = organize_files(reddir, table_path, 'astrometry',
                             excl_files, excl_objs, excl_filts)
    image_paths = file_df.paths
    
    # Prepare directories for processing
    mod_dir = maindir / 'processing' / 'mod-astro-input'  # Directory for modified images
    mod_dir.mkdir(parents=True, exist_ok=True)
    
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
            logger.info(f"Solution successful, saving WCS header to {output_path}")
            wcs_header.tofile(output_path, overwrite=True)
        else:
            logger.warning(f"Solution failed; skipping this image")

    # Ignore warnings about 'RADECSYS' header key being deprecated
    warnings.simplefilter("ignore", category=FITSFixedWarning)
    
    # Calibrate each image and collect output paths
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
    """
    Modify images to remove bad columns and create new FITS files with these changes.
    
    This function processes each image by masking out bad columns and other
    problematic regions, creating new modified FITS files in a specified directory.

    Parameters
    ----------
    image_paths : list of Path
        List of paths to the input images.
    mod_dir : Path
        Directory to save the modified FITS files.

    Returns
    -------
    list of Path
        List of paths to the modified FITS files.
    """
    logger.info("Zeroing out masked regions for faster astrometric solves")
    
    mod_paths = []
    for file in image_paths:
        mod_path = mod_dir / file.parent.name / (file.name.split('_')[0] + '_astroinput.fits')
        with fits.open(file) as hdul:
            data = np.array(hdul[0].data)
            try:
                # Create modified FITS files with masked regions set to 0
                mask = np.array(hdul['MASK'].data)
                mask = mask.astype(bool)
                data[mask] = 0
            except KeyError:
                # If no mask in FITS file, set bad columns to 0
                data[:, bad_columns] = 0
                logger.debug("No mask in FITS file--masking Nickel bad columns only")
            # Save the modified FITS file to the output directory
            fits.writeto(mod_path, data, hdul[0].header, overwrite=True)
        mod_paths.append(mod_path)
    
    return mod_paths


def get_astrometric_solves(image_paths, astro_dir, mode):
    """
    Retrieve local copies of astrometric solutions from previous runs.
    
    This function checks for previously solved images in the specified directory
    and returns their paths if they exist. It supports retrieving either the
    calibrated image or the source table with WCS information.

    Parameters
    ----------
    image_paths : list of Path
        List of paths to the input images.
    astro_dir : Path
        Path to the directory containing the astrometric solutions.
    mode : str
        Mode for retrieving files ('image' or 'corr').

    Returns
    -------
    list of Path
        List of paths to the local copies of astrometric solutions.
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
