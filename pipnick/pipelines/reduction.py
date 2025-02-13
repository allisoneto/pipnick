
from pathlib import Path
import warnings

from IPython import embed

import numpy as np

from astropy.nddata import CCDData
from astropy.stats import sigma_clip
from astropy.wcs.wcs import FITSFixedWarning
from astropy import units
import ccdproc

#from pipnick.utils.nickel_data import (gain, read_noise, bias_label, 
#                                       dome_flat_label, sky_flat_label,
#                                       sky_flat_label_alt,
#                                       dark_label, focus_label)

from pipnick import cameras
#from pipnick.utils.nickel_masks import get_masks_from_file
#from pipnick.utils.dir_nav import organize_files, norm_str
from pipnick.utils.dir_nav import build_metadata
from pipnick import logger


def reduce_all(rawdir, rdxdir=None, table=None, save=False, excl_files=None, excl_objs=None,
               excl_filts=None):
    """
    Perform reduction of raw astronomical data frames (overscan subtraction,
    bias subtraction, flat division, cosmic ray masking).

    Parameters
    ----------
    rawdir : str, Path
        Path to the parent directory of the raw directory
        containing the raw FITS files to be reduced.
    rdxdir : str, Path, optional
        Top-level directory for the reduced data.  If None, this is the parent
        directory of ``rawdir``.  I.e., if the raw directory is
        ``/User/janedoe/Nickel/raw``, ``rdxdir`` will be set to
        ``/User/janedoe/Nickel/``.
    table : str, Path, optional
        A file with the tabulated data to reduce.  If None, this file is
        automatically generated based on the data found and reduced in
        ``rawdir``; see :func:`~pipnick.utils.dir_nav.build_metadata`.
    save : bool, optional
        If True, save intermediate results during processing.
    excl_files : list, optional
        List of file stems to exclude (exact match not necessary).
    excl_objs : list, optional
        List of object strings to exclude (exact match not necessary).
    excl_filts : list, optional
        List of filter names to exclude.

    Returns
    -------
    list
        Paths to the reduced images written to disk.
    """
    # Setup paths
    _rawdir = Path(rawdir).absolute()
    if not _rawdir.is_dir():
        raise NotADirectoryError(f'{_rawdir} does not exist!')

    _rdxdir = _rawdir.parent if rdxdir is None else Path(rdxdir).absolute()
    if not _rdxdir.is_dir():
        _rdxdir.mkdir(parents=True)
    
    logger.info(f'Reducing data:')
    logger.info(f'    Raw data directory: {_rawdir}')
    logger.info(f'    Output directory: {_rdxdir}')

    # Organize raw files based on input directory or table
    metadata = build_metadata(_rawdir, rdx_table=table)
    nfiles = len(metadata)

    # Get the super bias
    super_bias = get_super_bias(metadata, save_dir=_rdxdir if save else None)

    # Get the super flat.  Try to get sky flats first:
    super_flat = get_super_flats(metadata, save_dir=_rdxdir if save else None, flattype='skyflat',
                                 bias=super_bias)

    # If there were no sky flats, try finding some dome flats
    if super_flat is None:
        super_flat = get_super_flats(metadata, save_dir=_rdxdir if save else None,
                                     flattype='domeflat', bias=super_bias)
    if super_flat is None:
        logger.warning('UNABLE TO BUILD FLAT-FIELDS!  Attempting to continue anyway...')
    
    # Assume all unknown frame types are on-sky science observations (but this
    # will include pointing and focus frames!)
    is_science = metadata['frametype'] == 'None'
    if np.sum(is_science) == 0:
        raise ValueError('No science frames found.')
    
    # Perform the basic processing of all science frames
    rdx_file = np.empty(len(metadata), dtype=object)
    for i in range(nfiles):
        if not is_science[i]:
            continue

        if super_flat is None or metadata['filter'][i] not in super_flat.keys():
            flat = None
        else:
            flat = super_flat[metadata['filter'][i]]

        # Trim and overscan correct it
        raw_file = Path(metadata['path'][i]).absolute() / metadata['file'][i]
        rdx_frame = process_frame(raw_file, flag_cosmics=True, bias=super_bias, flat=flat)
        rdx_file[i] = _rdxdir / f'{raw_file.stem}_rdx.fits'
        rdx_frame.write(rdx_file[i], overwrite=True)
        logger.info(f'Saved processed frame to {rdx_file[i]}')

    # Return the metadata table
    return metadata


def init_ccddata(frame, bpm=None, saturation=None, gain=None, readnoise=None):
    """
    Read data from a provided file.
    
    This function also uses the header information to identify the camera used
    to acquire the data.  If the camera is successfully identified, the
    camera information is used to 

        - apply any known bad-pixel mask
        - mask saturated pixels
        - convert the units from DN to counts by multiplying by the camera gain

    The relevant keyword arguments can be used to supplement information not
    available for a given camera or when the camera cannot be identified.  Note
    that values provided via a keyword *always* take precedence.

    Parameters
    ----------
    frame : str, Path
        Path to the FITS file.

    Returns
    -------
    CCDData
        Initialized and processed CCDData object.
    """
    warnings.simplefilter('ignore', FITSFixedWarning)
    _frame = Path(frame).absolute()

    logger.info(f'Initializing CCD imaging data in {_frame.name}.')
    ccd = CCDData.read(_frame, unit=units.adu)

    # Determine the camera used
    try:
        camera = cameras.identify_camera(ccd.header)
    except:
        warnings.warn(f'Unable to determine camera used for {_frame.name}.')
        camera = None
    else:
        ccd.meta['camera'] = camera
        camera = eval(f'cameras.{camera}')

    # Mask bad pixels
    _bpm = camera.bpm() if bpm is None and camera is not None else bpm
    if _bpm is not None:
        ccd.mask = _bpm

    # Mask saturated pixels; saturation is in ADU
    _sat = camera.saturation if saturation is None and camera is not None else saturation
    if _sat is not None:
        ccd.mask[ccd.data >= _sat] = True

    # Get and save the RN (in electrons) and gain (in electrons/ADU)
    ccd.meta['PROCGAIN'] = camera.gain if gain is None else gain
    ccd.meta['PROCRN'] = camera.readnoise if readnoise is None else readnoise

    # Apply gain if available
    if ccd.meta['PROCGAIN'] is not None:
        ccd.data = ccd.data * ccd.meta['PROCGAIN']
        ccd.unit *= (units.electron / units.adu)

    # Done
    return ccd


def process_frame(raw_frame, bias=None, flat=None, flag_cosmics=True):
    """
    Perform basic image processing

    This includes overscan subtraction, trimming, bias subtraction, and
    field-flattening.    

    Parameters
    ----------
    raw_frame : str, Path
        Data file to process.
    flag_cosmics : bool, optional
        Identify and mask cosmic rays
    bias : CCDData, optional
        Image to use for the bias subtraction
    flat : CCDData, optional
        Image to use for the field-flattening

    Returns
    -------
    CCDData
        Processed CCDData object with overscan subtracted and image trimmed.
    """
    warnings.simplefilter("error", RuntimeWarning)
    ccd = init_ccddata(raw_frame)
    try:
        oscansec = eval(f'cameras.{ccd.meta["camera"]}').oscansec(ccd.header)
    except:
        logger.warning('Could not determine overscan region!')
        oscansec = None
    else:
        logger.info('Subtracting overscan and trimming image.')
        ccd = ccdproc.subtract_overscan(ccd, fits_section=oscansec, overscan_axis=1)
        ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['DATASEC'])
    if bias is not None:
        logger.info('Subtracting super bias')
        ccd = ccdproc.subtract_bias(ccd, bias)
    if flat is not None:
        logger.info('Applying flat-field correction')
        ccd = ccdproc.flat_correct(ccd, flat)
    if flag_cosmics:
        logger.info(f'Removing cosmic rays from {raw_frame}')
        ccd = ccdproc.cosmicray_lacosmic(ccd, gain_apply=False, gain=1.0,
                                         readnoise=ccd.meta['PROCRN'], verbose=True)
    return ccd


def stack_frames(raw_frames, scale=False, bias=None, flat=None, flag_cosmics=True):
    """
    Process and stack a list of frames.

    Parameters
    ----------
    raw_frames : list
        List of file names to process and combine.
    scale : bool, optional
        Scale the image data by their (masked) mean before combining them?

    Returns
    -------
    CCDData
        Combined CCDData object.
    """
    # TODO: For *many* images, this may be too memory intensive.
    rdx_frames = [process_frame(frame, bias=bias, flat=flat, flag_cosmics=flag_cosmics) 
                    for frame in raw_frames]

    # TODO: Need to figure out how to handle the scaling; i.e., how is this
    # passed to the combiner?
    if scale:
        factor = np.array([float(f.mean().data) for f in rdx_frames])
        factor = factor[len(factor)//2] / factor
        rdx_frames = [f.multiply(s) for s,f in zip(factor, rdx_frames)]

    logger.info('Identifying outliers using sigma clipping')
    combiner = ccdproc.Combiner(rdx_frames)

    old_n_masked = 0
    new_n_masked = 1
    while new_n_masked > old_n_masked:
#        combiner.sigma_clipping(low_thresh=3, high_thresh=3, func=np.ma.mean)
        # TODO:
        #   - does ccdproc provide an iterative sigma clipper?
        #   - allow user to change the clipping parameters
        combiner.data_arr = sigma_clip(combiner.data_arr, sigma=5., maxiters=1,
                                       cenfunc=np.ma.mean, stdfunc=np.ma.std,
                                       axis=0, masked=True, copy=False)
        old_n_masked = new_n_masked
        new_n_masked = combiner.data_arr.mask.sum()

#    if frame_type == 'flat':
#        combiner.scaling = lambda arr: 1/np.ma.mean(arr)

    logger.info('Creating sigma-clipped average image.')
    mean = combiner.average_combine(scale_func=np.ma.mean)

    # TODO: Fill the masked regions?
    return mean


def get_super_bias(metadata, save_dir=None):
    """
    Create a combined (super) bias frame from individual bias frames.

    Parameters
    ----------
    metadata : Table
        Astropy Table containing file information.
    save_dir : Path or None, optional
        Directory to save the super bias frame.

    Returns
    -------
    CCDData
        SuperBias CCDData object.
    """
    # Select biases
    indx = metadata['frametype'] == 'bias'
    nbias = np.sum(indx)
    if nbias == 0:
        raise ValueError('No biases found in data files!')

    logger.info(f'Found {nbias} bias frames to process and combine.')

    frames = [Path(p) / f for p,f in zip(metadata['path'][indx], metadata['file'][indx])]
    super_bias = stack_frames(frames, flag_cosmics=False)

    if save_dir is None:
        return super_bias
    
    super_bias.header["OBJECT"] = "Bias"
    ofile = save_dir / 'Bias.fits'
    super_bias.write(ofile, overwrite=True)
    logger.info(f"Saved super bias to {ofile}")


def get_super_flats(metadata, save_dir=None, flattype='skyflat', bias=None, filter=None):
    """
    Create combined (super) flat frames (one per filter) from individual flat
    frames.

    Parameters
    ----------
    metadata : Table
        Astropy Table containing file information.
    save_dir : Path or None, optional
        Directory to save the super flat frames.
    flattype : str, optional
        Type of flats to use.  Must be 'domeflat' or 'skyflat'.
    bias : CCDData
        Super bias

    Returns
    -------
    dict
        Dictionary of super flat CCDData objects keyed by filter.
    """
    indx = metadata['frametype'] == flattype
    if np.sum(indx) == 0:
        logger.warning(f'No {flattype} frames found!')
        return None

    filters = np.unique(metadata['filter'][indx]) if filter is None else filter

    title = 'SkyFlat' if flattype == 'skyflat' else 'DomeFlat'

    super_flat = {}

    for f in filters:
        _indx = indx & (metadata['filter'] == f)
        if np.sum(_indx) == 0:
            continue
    
        logger.info(f'Processing and combining {np.sum(_indx)} {f} flats')
        frames = [Path(p) / f for p,f in zip(metadata['path'][_indx], metadata['file'][_indx])]
        super_flat[f] = stack_frames(frames, scale=True, bias=bias, flag_cosmics=False)
        if save_dir is None:
            continue

        super_flat[f].header["OBJECT"] = title
        ofile = save_dir / f'{title}_{f}.fits'
        super_flat[f].write(ofile, overwrite=True)
        logger.info(f'Saving {f} super {title} to {ofile}')

    return super_flat


def save_results(scifile_df, modifier_str, save_dir):
    """
    Save (partially) processed science files to the specified directory.

    Parameters
    ----------
    scifile_df : pd.DataFrame
        DataFrame containing processed science file information.
    modifier_str : str
        String to append to filenames to indicate processing stage.
    save_dir : Path
        Directory to save the processed files.

    Returns
    -------
    list
        List of paths to the saved files.
    """
    Path.mkdir(save_dir, exist_ok=True)
    logger.info(f"Saving {len(scifile_df.files)} _{modifier_str} images {save_dir.name} images to {save_dir}")
    save_paths = [save_dir / (path.stem.split('_')[0] + f"_{modifier_str}" + path.suffix) for path in scifile_df.paths]
    for file, path in zip(scifile_df.files, save_paths):
        file.write(path, overwrite=True)
    return save_paths


#bias_label = norm_str(bias_label)
#dome_flat_label = norm_str(dome_flat_label)
#sky_flat_label = norm_str(sky_flat_label)
#sky_flat_label_alt = norm_str(sky_flat_label_alt)
#dark_label = norm_str(dark_label)
#focus_label = norm_str(focus_label)

