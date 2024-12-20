
from pathlib import Path
import numpy as np
import warnings

from IPython import embed

from astropy.nddata import CCDData
from astropy.stats import sigma_clip
from astropy.wcs.wcs import FITSFixedWarning
from astropy import units
import ccdproc

from pipnick.utils.nickel_data import (gain, read_noise, bias_label, 
                                       dome_flat_label, sky_flat_label,
                                       sky_flat_label_alt,
                                       dark_label, focus_label)

from pipnick import cameras
#from pipnick.utils.nickel_masks import get_masks_from_file
from pipnick.utils.dir_nav import organize_files, norm_str
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
        ``rawdir``; see :func:`~pipnick.utils.dir_nav.organize_files`.
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
    logger.info(f"---- reduce_all() called on directory {rawdir}")
    _rawdir = Path(rawdir).absolute()
    if not _rawdir.is_dir():
        raise NotADirectoryError(f'{rawdir} does not exist!')

    _rdxdir = rawdir.parent if rdxdir is None else Path(rdxdir).absolute()
    if not _rdxdir.is_dir():
        _rdxdir.mkdir(parents=True)
    
    # Organize raw files based on input directory or table
    file_df = organize_files(_rawdir, table=table)
        #, 'reduction', excl_files, excl_objs, excl_filts)

    reddir = _rdxdir / 'reduced'
    procdir = _rdxdir / 'processing'
    reddir.mkdir(exist_ok=True)
    procdir.mkdir(exist_ok=True)
    
    # Initialize CCDData objects and remove cosmic rays
    logger.info("Initializing CCDData objects & removing cosmic rays")
    warnings.simplefilter("ignore", category=FITSFixedWarning)
#    ccd_objs = [init_ccddata(file) for file in file_dffiles]
#    file_df.files = ccd_objs
    
    # Get the super bias
    super_bias = get_super_bias(file_df, save=save, save_dir=procdir)

    # Get the super flat
    try:
        super_flat = get_super_flats(file_df, save=save, save_dir=procdir, flattype='skyflat',
                                    bias=super_bias)
    except:
        warnings.warn('Unable to construct sky flats.')
        super_flat = None
    # If there were no sky flats, try finding some dome flats
    if super_flat is None:
        try:
            super_flat = get_super_flats(file_df, save=save, save_dir=procdir, flattype='domeflat',
                                         bias=super_bias)
        except:
            warnings.warn('Unable to construct dome flats.')
    if super_flat is None:
        # TODO: Allow the code to keep going?
        raise ValueError('Cannot construct flats!')
    
    # Assume all unknown frame types are on-sky science observations (but this
    # will include pointing and focus frames!)
    indx = np.where(file_df['frametype'] == 'None')[0]
    if len(indx) == 0:
        raise ValueError('No science frames found.')

    reduced_science = []    
    for i in indx:
        # Trim and overscan correct it
        raw_frame = Path(file_df['path'][i]).absolute() / file_df['file'][i]
        rdx_frame = process_frame(raw_frame, flag_cosmics=True, bias=super_bias,
                                  flat=super_flat[file_df['filter'][i]])
        embed()
        exit()
        red_paths = save_results(raw_frame, rdx_frame, 'red')
        all_red_paths += red_paths        
        



    # Filter out non-science files
    scifiles_mask = ((file_df.objects != bias_label) &
                     (file_df.objects != dark_label) &
                     (file_df.objects != dome_flat_label) &
                     (file_df.objects != sky_flat_label) &
                     (file_df.objects != sky_flat_label_alt) &
                     (file_df.objects != focus_label)).values
    scifile_df = file_df.copy()[scifiles_mask]

    # Perform overscan subtraction & trimming
    logger.info(f"Performing overscan subtraction & trimming on {len(scifile_df.files)} science images")
    scifile_df.files = [trim_overscan(scifile) for scifile in scifile_df.files]
    if save_inters:
        save_results(scifile_df, 'over', procdir/'overscan')
    
    # Perform bias subtraction
    logger.info(f"Performing bias subtraction on {len(scifile_df.files)} science images")
    scifile_df.files = [ccdproc.subtract_bias(scifile, super_bias) 
                    for scifile in scifile_df.files]
    if save_inters:
        save_results(scifile_df, 'unbias', procdir/'unbias')

    # Perform flat division for each filter
    logger.info("Performing flat division")
    all_red_paths = []
    for filt in super_flats.keys():
        logger.debug(f"{filt} Filter:")
        scienceobjects = list(set(scifile_df.objects[scifile_df.filters == filt]))
        
        for scienceobject in scienceobjects:
            # Filter science files by object and filter
            sub_scifile_df = scifile_df.copy()[(scifile_df.objects == scienceobject) &
                                               (scifile_df.filters == filt)]
            # Create directory for each science target / filter combination
            sci_dir = reddir / (scienceobject + '_' + filt)
            
            # Perform flat division
            sub_scifile_df.files = [ccdproc.flat_correct(scifile, super_flats[filt]) 
                         for scifile in sub_scifile_df.files]
            
            red_paths = save_results(sub_scifile_df, 'red', sci_dir)
            all_red_paths += red_paths
    
    # Return
    logger.info(f"Fully reduced images saved to {reddir}")
    logger.info("---- reduce_all() call ended")
    return all_red_paths


def init_ccddata(frame, bpm=None, saturation=None, gain=None, readnoise=None):
    """
    Initialize a CCDData object from a FITS file and remove cosmic rays.

    Parameters
    ----------
    frame : str, Path
        Path to the FITS file.

    Returns
    -------
    CCDData
        Initialized and processed CCDData object.
    """
    _frame = Path(frame).absolute()

    logger.info(f'Loading data for {_frame.name}.')
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
    ccd.mask[ccd.data >= _sat] = True

    # Get and save the RN (in electrons) and gain (in electrons/ADU)
    ccd.meta['PROCGAIN'] = camera.gain if gain is None else gain
    ccd.meta['PROCRN'] = camera.readnoise if readnoise is None else readnoise

    # Apply gain
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
    oscansec = eval(f'cameras.{ccd.meta["camera"]}').oscansec(ccd.header)
    ccd = ccdproc.subtract_overscan(ccd, fits_section=oscansec, overscan_axis=1)
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['DATASEC'])
    if bias is not None:
        ccd = ccdproc.subtract_bias(ccd, bias)
    if flat is not None:
        ccd = ccdproc.flat_correct(ccd, flat)
    if flag_cosmics:
        logger.info(f'Removing cosmic rays from {raw_frame}')
        ccd = ccdproc.cosmicray_lacosmic(ccd, gain_apply=False, gain=1.0,
                                        readnoise=ccd.meta['PROCRN'], verbose=True)
    return ccd


def stack_frames(raw_frames, scale=False, bias=None, flat=None, flag_cosmics=True):
    """
    Stack frames by trimming overscan and combining them with sigma clipping.

    Parameters
    ----------
    raw_frames : list
        List of CCDData objects to combine.
    scale : bool, optional
        Scale the image data by their (masked) mean before combining them?

    Returns
    -------
    CCDData
        Combined CCDData object.
    """
    logger.info(f'Overscan subtracting and trimming {len(raw_frames)} images')
    rdx_frames = [process_frame(frame, bias=bias, flat=flat, flag_cosmics=flag_cosmics) 
                    for frame in raw_frames]
    
    if scale:
        factor = np.array([float(f.mean().data) for f in rdx_frames])
        factor = factor[len(factor)//2] / factor
        rdx_frames = [f.multiply(s) for s,f in zip(factor, rdx_frames)]

    logger.info('Sigma clipping images to combine')
    combiner = ccdproc.Combiner(rdx_frames)

    old_n_masked = 0
    new_n_masked = 1
    while new_n_masked > old_n_masked:
        combiner.data_arr = sigma_clip(combiner.data_arr, sigma=5., maxiters=1,
                                       cenfunc=np.ma.mean, stdfunc=np.ma.std,
                                       axis=0, masked=True, copy=False)
#        combiner.sigma_clipping(low_thresh=3, high_thresh=3, func=np.ma.mean)
        old_n_masked = new_n_masked
        new_n_masked = combiner.data_arr.mask.sum()

#    if frame_type == 'flat':
#        combiner.scaling = lambda arr: 1/np.ma.mean(arr)

    logger.info('Creating sigma-clipped average image.')
    mean = combiner.average_combine(scale_func=np.ma.mean)

    # TODO: Fill the masked regions?
    return mean


def get_super_bias(file_df, save=True, save_dir=None):
    """
    Create a combined (super) bias frame from individual bias frames.

    Parameters
    ----------
    file_df : Table
        DataFrame containing file information.
    save : bool, optional
        If True, save the super bias frame to disk.
    save_dir : Path or None, optional
        Directory to save the super bias frame.

    Returns
    -------
    CCDData
        SuperBias CCDData object.
    """
    # Select biases
    indx = file_df['frametype'] == 'bias'
    nbias = np.sum(indx)
    if nbias == 0:
        raise ValueError('No biases found in data files!')

    logger.info(f'Found {nbias} bias frames to combine.')

    frames = [Path(p) / f for p,f in zip(file_df['path'][indx], file_df['file'][indx])]
    super_bias = stack_frames(frames, flag_cosmics=False)
    if not save:
        return super_bias
    
    super_bias.header["OBJECT"] = "Bias"
    ofile = save_dir / 'Bias.fits'
    super_bias.write(ofile, overwrite=True)
    logger.info(f"Saved super bias to {ofile}")


def get_super_flats(file_df, save=True, save_dir=None, flattype='skyflat', bias=None):
    """
    Create combined (super) flat frames (one per filter) from individual flat
    frames.

    Parameters
    ----------
    file_df : Table
        DataFrame containing file information.
    save : bool, optional
        If True, save the super flat frames to disk.
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
    indx = file_df['frametype'] == flattype
    if np.sum(indx) == 0:
        raise ValueError(f'No {flattype} frames found!')
    filters = np.unique(file_df['filter'][indx])

    title = 'SkyFlat' if flattype == 'skyflat' else 'DomeFlat'

    super_flat = {}

    for f in filters:
        _indx = indx & (file_df['filter'] == f)
        if np.sum(_indx) == 0:
            continue
    
        logger.info(f'Processing and combining {np.sum(_indx)} {f} flats')
        frames = [Path(p) / f for p,f in zip(file_df['path'][_indx], file_df['file'][_indx])]
        super_flat[f] = stack_frames(frames, scale=True, bias=bias, flag_cosmics=False)
        if not save:
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


bias_label = norm_str(bias_label)
dome_flat_label = norm_str(dome_flat_label)
sky_flat_label = norm_str(sky_flat_label)
sky_flat_label_alt = norm_str(sky_flat_label_alt)
dark_label = norm_str(dark_label)
focus_label = norm_str(focus_label)

