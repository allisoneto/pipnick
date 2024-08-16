from pathlib import Path
import logging

from nickelpipeline.photometry.psf_photometry import psf_analysis, consolidate_groups
from nickelpipeline.photometry.aperture_photometry import aperture_analysis

logger = logging.getLogger(__name__)

def photometry_all(reddir, output_dir=None, thresh=8.0, group=False, mode='all',
                   fittype='circ', plot_final=False, plot_inters=False):
    r"""
    Perform photometric analysis on reduced files in a specified directory.

    Procedure:
    
        - Extract image files from the provided directory.
        - Perform PSF photometry on each image using 
          :func:`psf_analysis`.
        - Optionally consolidate detected sources into groups
          using :func:`consolidate_groups`.
        - Conduct aperture photometry on the results with
          :func:`aperture_analysis`.
        - Save the source catalogs to CSV files in the output directory.

    Parameters
    ----------
    reddir : :obj:`str`
        Path to directory containing reduced files for photometric analysis. 
        If the directory contains subdirectories, each subdirectory will be 
        processed independently.
    output_dir : :obj:`str`, optional
        Directory where the results will be saved. Defaults to a 
        ``/photometric/`` subdirectory at the same level as the parent 
        directory of ``reddir``.
    thresh : :obj:`float`, optional
        Threshold for source detection, defined as a multiple of the background 
        standard deviation. Default is 8.0.
    group : :obj:`bool`, optional
        If True, consolidate detected sources that belong to the same group 
        into a single source. Default is False.
    mode : :obj:`str`, optional
        Mode to run PSF photometry. Use 'all' for a complete analysis or 'new'
        for only analyzing newly detected sources. Choices are ['all', 'new']. 
        Default is 'all'.
    fittype : :obj:`str`, optional
        Type of Moffat fit to use for PSF photometry. Choices are 
        ['circ', 'ellip']. Default is 'circ'.
    plot_final : :obj:`bool`, optional
        If True, display final images with detected sources and their fluxes 
        labeled. Default is False.
    plot_inters : :obj:`bool`, optional
        If True, display intermediate images showing initial sources and 
        source groups for inspection. Default is False.

    Returns
    -------
    source_catalog_paths : `list`
        A list of file paths to the generated source catalogs (in CSV format).
        Each entry corresponds to an analyzed image file.
    """
    logger.info(f'---- photometry_all() called on directory {reddir}')
    
    # Convert reddir to a Path object and get all object directories
    reddir = Path(reddir)
    obj_dirs = [dir for dir in reddir.iterdir() if dir.is_dir()]
    # If no subdirectories are found, treat reddir as the only directory
    if len(obj_dirs) == 0:
        obj_dirs = [reddir]
    
    # Set default output directory if not provided
    if output_dir is None:
        output_dir = obj_dirs[0].parent.parent / 'photometric'
    # Define paths for different output categories
    unconsol_dir = output_dir / 'unconsolidated'
    consol_dir = output_dir / 'consolidated'
    proc_dir = output_dir / "proc_files"
    # Create necessary directories if they don't exist
    Path.mkdir(output_dir, exist_ok=True)
    Path.mkdir(unconsol_dir, exist_ok=True)
    Path.mkdir(consol_dir, exist_ok=True)
    Path.mkdir(proc_dir, exist_ok=True)
    
    source_catalog_paths = []  # List to hold paths to source catalogs

    # Process each object directory in reddir
    for obj_dir in obj_dirs:
        # Set the output directory based on whether grouping is enabled
        if group:
            output_dir = consol_dir / obj_dir.name
            output_file_append = '_photsrcs_consol.csv'
        else:
            output_dir = unconsol_dir / obj_dir.name
            output_file_append = '_photsrcs.csv'
        Path.mkdir(output_dir, exist_ok=True)
            
        # Process each file in the subdirectory
        for file in obj_dir.iterdir():
            # Set output file path
            filestem = file.stem.split('_')[0]
            output_file = output_dir / (filestem + output_file_append)
            
            # Perform PSF photometry on the file
            psf_data = psf_analysis(file, proc_dir, thresh=thresh, mode=mode, 
                                    fittype=fittype, plot_final=plot_final, 
                                    plot_inters=plot_inters)
            
            # Consolidate groups of sources if enabled
            if group:
                psf_data = consolidate_groups(psf_data)
            
            # Perform aperture photometry
            all_data = aperture_analysis(psf_data, file)
            
            # Save result
            all_data.write(output_file, format='csv', overwrite=True)
            source_catalog_paths.append(output_file)
    
    logger.info('---- photometry_all() call ended')
    # Return the list of paths to the source catalogs
    return source_catalog_paths