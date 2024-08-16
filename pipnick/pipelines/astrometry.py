import requests
import json
import time
import logging
from astropy.io import fits
from pathlib import Path
from pipnick.convenience.nickel_data import bad_columns
from pipnick.convenience.dir_nav import unzip_directories

######################################
####  ASTROMETRY.NET CALIBRATION  ####
######################################

logger = logging.getLogger(__name__)

def astrometry_all(reddir, api_key, output_dir=None, mode='image', resolve=False):
    """
    Runs astrometry.net calibration on all images input, and saves the calibrated
    fits files to output_dir. Uses astrometry.net's API.

    Args:
        reddir (Path or str): path to directory containing images to calibrate
        output_dir (str): path to output image folder
        mode (str): Whether to return paths to calibrated image or source table. 'image' or 'corr'.
        resolve (bool): If True, re-solves images with previously generated local solves

    Returns:
        list: list of relative paths (str) to all calibrated fits images
    """
    logger.info(f'---- astrometry_all() called on directory {reddir}')
    reddir = Path(reddir)
    red_files = unzip_directories([reddir], output_format='Path')
    if output_dir is None:
        output_dir = red_files[0].parent.parent.parent / 'astrometric'
    
    logger.info('---- astrometry_all() call ended')
    return run_astrometry(red_files, api_key, output_dir, mode, resolve)
    
    
def run_astrometry(image_paths, api_key, output_dir, mode='image', resolve=False):
    """
    Runs astrometry.net calibration on all images input, and saves the calibrated
    fits files to output_dir. Uses astrometry.net's API.

    Args:
        image_paths (list): list of paths (str) to all images
        output_dir (str): path to output image folder
        mode (str): Whether to return paths to calibrated image or source table. 'image' or 'corr'.
        resolve (bool): If True, re-solves images with previously generated local solves

    Returns:
        list: list of relative paths (str) to all calibrated fits images
    """
    if not isinstance(image_paths, list):
        image_paths = [Path(image_paths)]
    else:
        image_paths = [Path(path) for path in image_paths]
    
    # Makes output folder if it doesn't already exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    newimg_dir = output_dir / 'astroimg'
    newimg_dir.mkdir(parents=True, exist_ok=True)
    corr_dir = output_dir / 'corr'
    corr_dir.mkdir(parents=True, exist_ok=True)
    
    # Modify images to remove the Nickel Telescope's bad columns
    logger.info("Zeroing out masked regions for faster astrometric solves")
    mod_dir = output_dir.parent / 'raw-astro-input'  # Directory for modified images
    mod_dir.mkdir(parents=True, exist_ok=True)
    
    mod_paths = []
    for file in image_paths:
        mod_path = mod_dir / file.name
        if not mod_path.exists():
            with fits.open(file) as hdul:
                data = hdul[0].data
                try:
                    # Creating modified FITS files w/ masked regions set to 0
                    mask = hdul['MASK'].data
                    data[mask] = 0
                except KeyError:
                    # If no mask in FITS file, sets bad columns = 0
                    data[:, bad_columns] = 0
                    logger.debug("No mask in FITS file--masking Nickel bad columns")
                # Save the modified FITS file to the output directory
                hdul.writeto(mod_path, overwrite=True)
        mod_paths.append(mod_path)
    
    # Log in to Nova and get session key
    logger.info("Connecting to astrometry.net")
    R = requests.post(
        "http://nova.astrometry.net/api/login",
        data={"request-json": json.dumps({"apikey": str(api_key)})},  # API key can be found in your account -- currently using Allison's
    )
    dictionary = json.loads(R.text)
    session_key = dictionary["session"]
    
    # Wrapper function to convert a single photo (to be used on all photos below)
    def convert(image_path):
        logger.info(f"Submitting image {image_path.name} to astrometry.net")
        # Use session key to upload file(s) to Astrometry
        url = "http://nova.astrometry.net/api/upload"
        files = {"file": (image_path.name, open(image_path, "rb"))}
        data = {
            "request-json": '{"publicly_visible": "y", "allow_modifications": "d", "session": "' + session_key + '", "allow_commercial_use": "d"}' #, "scale_type":"ul", "scale_lower":0.05, "scale_upper":5}'
        }
        R = requests.post(url, files=files, data=data)
        dictionary = json.loads(R.text)
        subid = dictionary['subid']

        # Get jobid from submission status
        # Requires checking every 5 sec until requests.post actually has a JobID
        # Otherwise an error will be thrown on some images
        logger.info("Waiting for astrometry.net to accept image. May take several minutes")
        url = f"http://nova.astrometry.net/api/submissions/{subid}"
        jobid_found = False
        while jobid_found is False:
            logger.debug("Waiting for JobID")
            time.sleep(10)
            R = requests.post(url)
            dictionary = json.loads(R.text)
            logger.debug(dictionary)
            try:
                jobid = dictionary["jobs"][0]
                if jobid is None:
                    continue
                jobid_found = True
            except:
                continue
        
        logger.info("Image accepted--beginning solve.")
        logger.debug(f"JobID = {jobid}")

        # Every t (currently 15) seconds, check if image has been calibrated
        # Exit & fail if taking > 10 minutes to calibrate
        logger.info(f"Time elapsed = {0} seconds")
        url = f"http://nova.astrometry.net/api/jobs/{jobid}"
        def time_check(t, interval):
            time.sleep(interval)
            t += interval
            logger.info(f"Time elapsed = {t} seconds")
            R = requests.post(url)
            data = json.loads(R.text)
            logger.debug(data)
            if data['status'] == 'success':
                logger.info('Job succeeded, getting calibrated fits file')
                return True
            elif data['status'] == 'failed':
                logger.warning('Job failed')
                return False
            if t > 225:
                logger.warning('Maximum time elapsed... exiting')
                return False
            return time_check(t, interval)
        
        success = time_check(0, 15)
        logger.debug(f"success status = {success}")
        
        if success:
            logger.debug(f"Trying to save calibrated image & corr")
            # Create path to image folder
            output_path_stem = image_path.stem.split('_')[0]
            
            # Get FITS file & save to output path (in folder just created ^)
            url_image = f"http://nova.astrometry.net/new_fits_file/{jobid}"
            url_corr = f"http://nova.astrometry.net/corr_file/{jobid}"
            output_path_img = newimg_dir / f"{output_path_stem}_astro.fits"
            output_path_corr = corr_dir / f"{output_path_stem}_astro_corr.fits"
            
            new_image = requests.get(url_image)
            corr_table = requests.get(url_corr)
            
            with open(output_path_img, 'wb') as img:
                img.write(new_image.content)
            with open(output_path_corr, 'wb') as corr:
                corr.write(corr_table.content)
            
            # Return output_path for other functions to find this file
            logger.info(f"Calibrated image & corr saved to {output_path_img} (corr w/ _corr.fits)")
            if mode == 'image':
                return output_path_img
            elif mode == 'corr':
                return output_path_corr
        else:
            logger.debug("Saving nothing")
            return None
    
    # Calibrate each image piece and collect output_paths
    calibrated_fits_paths = []
    for image_path in mod_paths:
        output_path_stem = image_path.stem.split('_')[0]
        if mode == 'image':
            output_path = newimg_dir / f"{output_path_stem}_astro.fits"
        elif mode == 'corr':
            output_path = corr_dir / f"{output_path_stem}_astro_corr.fits"
        if not resolve and output_path.exists():
            logger.info(f"Returning local copy of {image_path.name}'s solution; astrometry.net not used")
            calibrated_fits_paths.append(output_path)
        else:
            try:
                calib_fits = convert(image_path)
                if calib_fits is not None:
                    calibrated_fits_paths.append(calib_fits)
            except requests.exceptions.ConnectionError:
                logger.warning(f"***Connection Error encountered; skipping image {image_path.name} & waiting 15 sec***")
                time.sleep(15)
        # os.remove(image_path)
    return calibrated_fits_paths


def get_astrometric_solves(image_paths, output_dir, mode):
    """
    Returns any local copies of astrometric solves stored from previous runs of
    run_astrometry(). Skips if image has not yet been solved.

    Args:
        image_paths (list): list of relative paths (str) to all images
        output_dir (str): path to output image folder
        mode (str): Whether to return paths to calibrated image or .corr file w/ source table

    Returns:
        list: list of relative paths (str) to all calibrated fits images
    """
    
    logger.info("Returning local copies of astrometric solves; astrometry.net not used")
    calibrated_fits_paths = []
    output_dir = Path(output_dir)
    for image_path in image_paths:
        output_path_stem = Path(image_path).stem
        if mode == 'image':
            output_path = output_dir / f"{output_path_stem[:-5]}.fits"
        elif mode == 'corr':
            output_path = output_dir / f"{output_path_stem[:-5]}_corr.fits"
        if output_path.exists():
            logger.debug(f"Found calibrated image {Path(image_path).name}; appending to list")
            calibrated_fits_paths.append(output_path)
    return calibrated_fits_paths
