import numpy as np
from pathlib import Path
import logging

from nickelpipeline.psf_analysis.moffat.stamps import generate_stamps
from nickelpipeline.psf_analysis.moffat.fit_psf import fit_psf_single, fit_psf_stack, psf_plot
from nickelpipeline.convenience.dir_nav import unzip_directories
from nickelpipeline.psf_analysis.moffat.model_psf import FitEllipticalMoffat2D, FitMoffat2D
from nickelpipeline.convenience.nickel_data import plate_scale_approx


logger = logging.getLogger(__name__)


def get_source_pars(path_list, category_str=None, fittype='ellip'):
    """
    Extract source coordinates and Moffat fit parameters from image data.

    Parameters
    ----------
    path_list : list
        List of paths (directories or files) to unzip.
    category_str : str, optional
        Category string for identifying the path to data.
    fittype : str, optional
        Type of model to fit ('ellip' or 'circ').

    Returns
    -------
    source_coords : `numpy.ndarray`
        Array of source coordinates.
    source_pars : `numpy.ndarray`
        Array of source parameters (x0, y0, amplitude, gamma1, gamma2, phi, alpha, background).
    img_nums : `numpy.ndarray`
        Array of image numbers corresponding to the sources.
    """
    # Unzip directories to get image files
    images = unzip_directories(path_list, output_format='Path')

    # Create output directories
    proc_dir = Path('.').resolve() / "proc_files"
    Path.mkdir(proc_dir, exist_ok=True)
    proc_subdir = proc_dir / fittype
    Path.mkdir(proc_subdir, exist_ok=True)
    base_parent = proc_subdir / category_str
    Path.mkdir(base_parent, exist_ok=True)
    base = proc_subdir / category_str / category_str

    # Generate stamps (image of sources) for image data
    generate_stamps(images, output_base=base)

    # Fit PSF models and get source coordinates and parameters
    source_coords, source_fits, img_nums = fit_psf_single(base, len(images))
    source_pars = np.array([fit.par for fit in source_fits])
    return source_coords, source_pars, img_nums


def get_graphable_pars(file_paths, group_name, verbose=False):
    """
    Fit PSF and extract parameters for given files (all stars in these files are
    stacked), storing intermediates in folder proc_files/elliptical/group_name.

    Parameters
    ----------
    file_paths : list
        List of file paths to analyze.
    group_name : str
        Folder for storing intermediates (proc_files/group_name).
    verbose : bool, optional
        If True, print detailed output during processing.

    Returns
    -------
    fwhm : float
        FWHM.
    ecc : float
        FWHM eccentricity.
    phi : float
        Rotation angle phi.
    """
    # Define directory and base path for processed files
    proc_dir = Path('.').resolve() / "proc_files"
    Path.mkdir(proc_dir, exist_ok=True)
    proc_subdir = proc_dir / 'ellip'
    Path.mkdir(proc_subdir, exist_ok=True)
    base_parent = proc_subdir / group_name
    Path.mkdir(base_parent, exist_ok=True)
    base = proc_subdir / group_name / group_name

    # Generate image stamps for the given files
    generate_stamps(file_paths, output_base=base)

    # Fit PSF stack and get the fit results
    psf_file = Path(f'{str(base)}.psf.fits').resolve()  # PSF info stored here
    fit = fit_psf_stack(base, 1, fittype='ellip', ofile=psf_file)

    # Plot PSF and get FWHM and phi values
    plot_file = Path(f'{str(base)}.psf.pdf').resolve()  # Plots stored here
    psf_plot(plot_file, fit, fittype='ellip', plot_fit=True)
    fwhm = get_param_list('fwhm', np.array([fit.par]), (1,))[0][0]
    ecc = get_param_list('ecc', np.array([fit.par]), (1,))[0][0]
    phi = get_param_list('phi', np.array([fit.par]), (1,))[0][0]

    if verbose:
        print(f"Avg FWHM = {fwhm:3f}")
        print(f"FWHM_ecc = {ecc:3f}")
        print(f"Rotation angle phi = {phi:3f}")

    return fwhm, ecc, phi


def get_param_list(param_type, pars, shape, img_nums=None):
    """
    Generate the desired single parameter list, color range, and title
    for contour plotting based on Moffat pars.

    Parameters
    ----------
    param_type : str
        Type of parameter ('fwhm', 'fwhm residuals', 'phi', 'ecc').
    pars : `numpy.ndarray`
        Fit parameters (list of par).
    shape : tuple
        Shape to output param_list.
    img_nums : `numpy.ndarray`, optional
        Image number for each source.

    Returns
    -------
    param_list : `numpy.ndarray`
        List of parameter values.
    color_range : list
        Range of colors for plotting.
    title : str
        Title for the plot.

    Raises
    ------
    ValueError
        If the input `param_type` is not 'fwhm', 'phi', 'ecc', or 'fwhm residuals'.
    """
    if param_type == 'fwhm':
        # Calculate FWHM (average between semi-major and minor axes)
        param_list = (FitMoffat2D.to_fwhm(pars[:, 3], pars[:, 6]) +
                      FitMoffat2D.to_fwhm(pars[:, 4], pars[:, 6])) / 2 * plate_scale_approx
        color_range = [1.5, 2.7]  # Optimized for Nickel 06-26-24 data
        title = "FWHM (arcsec)"
    elif param_type == 'fwhm residuals':
        # Calculate FWHM residual (relative to minimum FWHM in image)
        fwhm_list = (FitMoffat2D.to_fwhm(pars[:, 3], pars[:, 6]) +
                     FitMoffat2D.to_fwhm(pars[:, 4], pars[:, 6])) / 2
        mins = {img_num: np.min(fwhm_list[img_nums == img_num])
                for img_num in list(set(img_nums))}
        param_list = np.array([fwhm_list[i] - mins[img_num]
                               for i, img_num in enumerate(img_nums)]) * plate_scale_approx
        color_range = [0.0, 0.36]
        title = "FWHM Residuals (arcsec)"
    elif param_type == 'phi':
        # Convert phi rotation angle relative to x-axis from the original phi
        param_list = np.array([FitEllipticalMoffat2D.get_nice_phi(smooth_par)
                               for smooth_par in pars])
        color_range = [-45., 45.]
        title = "Phi Rotation Angle (deg)"
    elif param_type == 'ecc':
        # Calculate eccentricity
        param_list = []
        for smooth_par in pars:
            fwhm1 = FitMoffat2D.to_fwhm(smooth_par[3], smooth_par[6])
            fwhm2 = FitMoffat2D.to_fwhm(smooth_par[4], smooth_par[6])
            param_list.append(np.sqrt(np.abs(fwhm1 ** 2 - fwhm2 ** 2)) / max(fwhm1, fwhm2))
        param_list = np.array(param_list)
        color_range = [0.29, 0.65]  # Optimized for Nickel 06-26-24 data
        title = "Eccentricity"
    else:
        raise ValueError("Input param_type must be 'fwhm' or 'phi'")

    param_list = param_list.reshape(shape)
    return param_list, color_range, title