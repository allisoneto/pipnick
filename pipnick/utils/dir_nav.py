import logging
from pathlib import Path
from pipnick.utils.fits_class import Fits_Simple

logger = logging.getLogger(__name__)


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