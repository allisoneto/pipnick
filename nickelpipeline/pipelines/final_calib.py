import logging
from pathlib import Path
from nickelpipeline.final_calibration.combo_astro_photo import convert_coords_all

#################################################
####  This pipeline is currently incomplete  ####
#################################################

logger = logging.getLogger(__name__)

def final_calib_all(photo_dir, astro_dir, output_dir=None):
    """
    Perform the final calibration step by converting pixel coordinates in
    photometric source catalogs to include sky coordinates (RA/Dec).

    This function takes directories containing photometric source catalogs and 
    astrometric calibration files, processes them to generate calibrated catalogs 
    with RA/Dec coordinates using astrometric calibration data, and saves the
    results to an output directory.

    Parameters
    ----------
    photo_dir : Path or str
        Directory containing photometric source catalogs.
    astro_dir : Path or str
        Directory containing astrometric calibration files.
    output_dir : Path or str, optional
        Directory where the final calibrated catalogs will be saved. 
        If not provided, the output directory will be set to a default location 
        relative to `astro_dir`.

    Returns
    -------
    dict
        A dictionary where keys are object names and values are the corresponding 
        calibrated catalogs with RA/Dec coordinates.
    """
    
    photo_dir = Path(photo_dir)
    astro_dir = Path(astro_dir)
    
    # Set default output directory if none is provided
    if output_dir is None:
        ouptut_dir = astro_dir.parent.parent / 'final_calib'
    Path.mkdir(ouptut_dir, exist_ok=True)
    
    # Convert coordinates for all photometric catalogs and save the results
    astrophot_datas = convert_coords_all(photo_dir, astro_dir, ouptut_dir)
    
    return astrophot_datas
