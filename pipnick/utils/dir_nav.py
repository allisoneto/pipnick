import logging
from pathlib import Path
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from pipnick.utils.fits_class import Fits_Simple
from pipnick.utils.nickel_data import focus_label


logger = logging.getLogger(__name__)


def organize_files(datadir, table_path, mode,
                   excl_files, excl_objs, excl_filts):
    """
    Extract, organize files by metadata, and apply exclusions to
    produce a pandas DataFrame of images to perform functions like
    reduction, astrometry, photometry, and final_calibration on.
    Saves information to or draws information from a table file, and 
    comments out files to be excluded.

    Parameters
    ----------
    datadir : str or Path
        Path to the directory containing the FITS files to be analyzed.
    table_path : str or Path, optional
        Path to a pipnick-specific table file with information about
        which raw FITS files to process. Must be produced by organize_files()
    mode : str
        Function for which organize_files() is run, & name of new table output
        ('reduction', 'astrometry', 'photometry', 'final_calibration')
    excl_files : list
        List of file stems to exclude (exact match not necessary).
    excl_objs : list
        List of object strings to exclude (exact match not necessary).
    excl_filts : list
        List of filter names to exclude.

    Returns
    -------
    pd.DataFrame
        DataFrame containing organized file information.
    """
    
    if table_path is None:
        # Set path to save table of files processed in reduction
        table_path = datadir.parent / f'{mode}_files.tbl'
        
        # Extract files from the specified directory
        if mode == 'reduction':
            files = [file for file in datadir.iterdir() if file.is_file()]
        else:
            obj_dirs = [dir for dir in datadir.iterdir() if dir.is_dir()]
            files = [file for obj_dir in obj_dirs for file in obj_dir.iterdir()]
        logger.info(f"{len(files)} files extracted from directory {datadir}")
        
        # Create DataFrame with file metadata
        obj_list = []
        filt_list = []
        for file in files:
            hdul = fits.open(str(file))
            obj_list.append(norm_str(hdul[0].header["OBJECT"]))
            filt_list.append(hdul[0].header["FILTNAM"])
            hdul.close()

        file_df = pd.DataFrame({
            "names": [file.stem for file in files],
            "files": files,
            "objects": obj_list,
            "filters": filt_list,
            "paths": files
            })
    
        # Save the table file for future reference
        logger.info(f"Saving table of {mode} file data to {table_path}")
        logger.info(f"You can change the file name to {mode}_files.yml for ease of commenting out files in VS Code (ctrl + '/')")
        file_table = Table.from_pandas(file_df)
        file_table.remove_column('files')
        file_table.write(table_path, format='ascii.fixed_width', overwrite=True)
    else:
        # Extract files from an astropy Table file
        logger.info(f"Files will be extracted from Astropy table file {table_path}, not directory {datadir}")
        file_table = Table.read(table_path, format='ascii.fixed_width')
        # Convert astropy table to pandas DataFrame
        file_df = file_table.to_pandas()
        file_df.insert(1, "files", file_df.paths)
        file_df.paths = [Path(file_path) for file_path in file_df.paths]
        logger.info(f"{len(file_df.paths)} files extracted from table file")
    
    # Apply manual exclusions based on provided criteria
    # (excludes file if any str in excl_list is in the the file's excl_type)
    excl_types = ['name', 'object', 'filter', 'object']
    excl_lists = [excl_files, excl_objs, excl_filts, [focus_label]]
    excl_file_names_all = []
    for excl_type, excl_list in zip(excl_types, excl_lists):
        excl_func = create_exclusion_func(excl_list)
        axis = file_df[excl_type+'s']
        excl_file_names = list(file_df.names[axis.apply(lambda x: not excl_func(x))])
        excl_file_names_all += excl_file_names
        file_df = file_df[axis.apply(excl_func)]
        if len(excl_file_names) > 0:
            logger.info(f"Manually excluding files with {excl_list} in {excl_type}: {excl_file_names}")
    
    # Add '#' to excluded files' rows in table file to ignore them in future
    already_excl_lines = comment_out_rows(excl_file_names_all, table_path, modify=True)
    if len(already_excl_lines) > 0:
        logger.info(f"Automatically excluding files already commented out in the table file: {already_excl_lines}")
    if len(excl_file_names_all) > 0:
        logger.info(f"Modifying table at {table_path} to ignore manually excluded files {excl_file_names_all} in future")
    
    # Return
    return file_df


def comment_out_rows(excluded_file_names, table_file, modify=True):
    """
    Comment out specified rows in a table file based on exclusion criteria.

    Parameters
    ----------
    excluded_file_names : list
        List of file names to comment out.
    table_file : str
        Path to the table file to modify.
    modify : bool, optional
        Whether to modify the file or not.

    Returns
    -------
    list
        List of file names that were already commented out.
    """
    with open(table_file, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    already_excl_lines = []
    for line in lines:
        if line.strip().startswith('#'):
            # Retrieves just the file stem (i.e. 'd1001')
            already_excl_lines.append(line.split('|')[1].split(' ')[1])
            new_lines.append(line)
        elif modify and any(file_name in line for file_name in excluded_file_names):
            new_lines.append('#' + line)    # "comments out" the row
        elif modify:
            new_lines.append(line)
    
    if modify:
        with open(table_file, 'w') as f:
            f.writelines(new_lines)
    
    return already_excl_lines


def norm_str(s):
    """
    Normalize a string for comparison purposes--all caps, no spaces.
    'Sky flat' -> 'SKYFLAT'

    Parameters
    ----------
    s : str or list
        String or list of strings to normalize.

    Returns
    -------
    str or list
        Normalized string or list of normalized strings.
    """
    if isinstance(s, list):
        return [norm_str(elem) for elem in s]
    return s.upper().replace(' ', '')

def create_exclusion_func(exclude_list):
    """
    Create a function to determine if a file should be excluded based on a list of criteria.

    Parameters
    ----------
    exclude_list : list
        List of criteria for exclusion.

    Returns
    -------
    function
        Function that takes a target (string) and returns True if it should be excluded.
    """
    if exclude_list is None:
        return lambda _: True
    exclude_list = [norm_str(obj_str) for obj_str in exclude_list]
    def excl_func(target):
        target_str = norm_str(target)
        is_excluded = any(excluded_str in target_str for excluded_str in exclude_list)
        return not is_excluded
    return excl_func


def unzip_directories(path_list: list, output_format: str = 'Fits_Simple', allow_exceptions: bool = False) -> list:
    """
    Extract all files in a list of directories or files (or their subdirectories) and
    return a list of image objects.

    Parameters
    ----------
    path_list : list
        A list of paths to directories or files that need to be unzipped.
    output_format : str, optional
        The format of the output objects. Can be 'Fits_Simple' (default) or 'Path'.
    allow_exceptions : bool, optional
        Whether to allow and handle exceptions for file-related errors such
        as `KeyError` or `OSError`. Defaults to False.

    Returns
    -------
    list
        A list of image objects in the specified output format. The objects can be
        instances of `Fits_Simple` or `Path`, depending on input `output_format`.
    """
    if output_format == 'Path':
        output = Path
    elif output_format == 'Fits_Simple':
        output = Fits_Simple
    
    images = []
    for elem_path in path_list:
        elem_path = Path(elem_path)
        if elem_path.is_dir():
            for sub_elem_path in elem_path.iterdir():
                images += unzip_directories([sub_elem_path], output_format, allow_exceptions)
        elif elem_path.is_file():
            if allow_exceptions:
                try:
                    images.append(output(elem_path))
                except (KeyError, OSError):
                    logger.debug(f"unzip_directories() threw KeyError or OSError on {elem_path}")
            else:
                images.append(output(elem_path))
    return images