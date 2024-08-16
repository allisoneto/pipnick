
#################################################################
########  Locate sources and create "stamps" to analyze  ########
#################################################################

import numpy as np
from pathlib import Path
from matplotlib import pyplot, ticker
from matplotlib.backends.backend_pdf import PdfPages
import logging

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from skimage import transform

from nickelpipeline.convenience.fits_class import Fits_Simple

logger = logging.getLogger(__name__)


def show_stamps(stamps, pdf_file, vmin=None, vmax=None):
    """
    Display and save the stamps in a PDF file.

    Args:
        stamps (ndarray): Array of stamp images.
        pdf_file (str): Path to the output PDF file.
        vmin (float, optional): Minimum data value to map to colormap.
        vmax (float, optional): Maximum data value to map to colormap.
    """
    dim = 8  # Dimension of the grid for displaying stamps
    start = 0.01  # Starting position for the grid
    buff = 0.01  # Buffer between stamps
    dw = (1-2*start-(dim-1)*buff)/dim  # Width of each stamp
    nstamps = stamps.shape[0]  # Total number of stamps

    # Create a PDF file to save the stamps
    with PdfPages(pdf_file) as pdf:
        i = 0
        while i < nstamps:
            w, h = pyplot.figaspect(1)
            fig = pyplot.figure(figsize=(1.5*w, 1.5*h))
            for j in range(dim):
                for k in range(dim):
                    if i < nstamps:
                        # Add a subplot for each stamp
                        ax = fig.add_axes([start+k*(dw+buff), start+(dim-j-1)*dw+(dim-j-1)*buff,
                                           dw, dw])
                        ax.xaxis.set_major_formatter(ticker.NullFormatter())
                        ax.xaxis.set_visible(False)
                        ax.yaxis.set_major_formatter(ticker.NullFormatter())
                        ax.yaxis.set_visible(False)
                        ax.imshow(stamps[i], origin='lower', interpolation='nearest',
                                  vmin=vmin, vmax=vmax)
                        ax.text(0.1, 0.9, str(i), ha='center', va='center', 
                                transform=ax.transAxes, color='C1')
                        ax.scatter([stamps.shape[1]//2], [stamps.shape[2]//2], 
                                   marker='+', s=50, color='C3', zorder=10, lw=0.5)
                        i += 1
            pdf.savefig()  # Save the figure to the PDF file
            fig.clear()
            pyplot.close()

def generate_stamps(images, output_base, thresh=15.):
    """
    Generate and analyze stamps from the given images.

    Args:
        images (list): List of image file paths.
        output_base (Path): Base of path to files for storing stamp data
        thresh (float): Star detection threshold (threshold = thresh*std)
    """
    num_images = len(images)
    
    # Load images and create masked arrays
    images = [Fits_Simple(image) for image in images]
    masked_images = np.ma.array([image.masked_array for image in images])
    
    ofits = output_base.with_suffix('.rdx.fits')
    src_ofile = output_base.with_suffix('.src.db')

    # Setting up stamp & source detection parameters
    default_fwhm = 7.0
    stamp_width = int(np.ceil(5 * default_fwhm))
    if stamp_width % 2 == 0:
        stamp_width += 1
    stamp_shape = (stamp_width, stamp_width)
    x, y = np.meshgrid(np.arange(stamp_width, dtype=float) - stamp_width // 2, 
                       np.arange(stamp_width, dtype=float) - stamp_width // 2)
    r = np.sqrt(x**2 + y**2)

    mean_rkron = np.zeros(num_images, dtype=float)  # Mean Kron radius for each image
    source_data = np.zeros((0, 9), dtype=float)  # Array to hold source data
    all_stamps = np.empty(num_images, dtype=object)  # Array to hold all stamps

    # Detect sources and generate stamps
    for index in range(num_images):
        # Source detection
        mean, median, std = sigma_clipped_stats(masked_images[index], sigma=3.)
        starfind = DAOStarFinder(fwhm=default_fwhm, threshold=thresh*std)
        sources = starfind(masked_images[index].filled(0.0) - median, mask=np.ma.getmaskarray(masked_images[index]))
        if sources is None:
            continue
        nsources = len(sources)
        ap_centers = np.column_stack((sources['xcentroid'], sources['ycentroid']))  # Source positions

        # Stamp setup
            # Scikit-Image uses the bottom left corner to define the coordinate system.  
            # Photutils uses the center of the pixel. ...  I think.  In any case, there 
            # seems to be an 0.5 pixel offset, which is why the transformation below 
            # has that additional offset
        center = np.atleast_1d(stamp_shape) / 2  # Center of the stamp
        stamps = np.ma.zeros((nsources,) + stamp_shape, dtype=float)  # Array to hold the stamps
        rkron = np.zeros(nsources, dtype=float)  # Array to hold Kron radii
        kflux_raw = np.zeros(nsources, dtype=float)  # Array to hold raw flux
        bkg = np.zeros(nsources, dtype=float)  # Array to hold background values

        # Process each detected source and add to stamp
        for i in range(nsources):
            t = transform.EuclideanTransform(translation=-ap_centers[i][::-1] - 0.5 + center)
            stamps[i] = transform.warp(masked_images[index].filled(0.0).T, t.inverse, output_shape=stamp_shape, cval=0, order=0).T
            stamps[i, stamps[i] == 0.] = np.ma.masked  # Mask the background
            if np.all(np.ma.getmaskarray(stamps[i, r > 3 * default_fwhm])):
                continue
            bkg[i] = np.ma.median(stamps[i, r > 3 * default_fwhm])
            stamps[i] -= bkg[i]  # Subtract the background
            if np.all(np.ma.getmaskarray(stamps[i])):
                continue
            _rkron = np.ma.sum(r * stamps[i]) / np.ma.sum(stamps[i])
            if _rkron == np.ma.masked or _rkron < 0:
                continue
            rkron[i] = _rkron  # Calculate the Kron radius
            indx = r < rkron[i]
            if not np.any(indx):
                continue
            kflux_raw[i] = np.ma.sum(stamps[i, indx])  # Calculate the raw flux

        mean_rkron[index] = sigma_clipped_stats(rkron, sigma=3., maxiters=None)[0]
        kflux = np.array([np.sum(s[r < mean_rkron[index]]) for s in stamps])

        # Remove sources that are too close to each other
        i, j = np.triu_indices(nsources, k=1)
        dist = np.sqrt((sources['xcentroid'][j] - sources['xcentroid'][i])**2 + 
                       (sources['ycentroid'][j] - sources['ycentroid'][i])**2)
        too_close = np.identity(nsources).astype(bool)
        min_dist_fwhm = 5
        too_close[i, j] = dist < min_dist_fwhm * default_fwhm
        too_close[j, i] = dist < min_dist_fwhm * default_fwhm

        # Remove sources with low flux
        min_flux = 1e3
        _kflux = np.ma.MaskedArray(np.tile(kflux, (nsources, 1)), mask=np.logical_not(too_close))
        keep = ((np.ma.argmax(_kflux, axis=1) == np.arange(nsources)) & 
                (kflux > min_flux) & 
                (np.sum(stamps.mask, axis=(1, 2)) == 0))
        nkeep = np.sum(keep)
        if nkeep == 0:
            all_stamps[index] = None
            continue

        srt = np.argsort(kflux[keep])[::-1]
        all_stamps[index] = stamps[keep][srt]
        source_data = np.append(source_data,
                                np.column_stack((np.full(nkeep, index+1, dtype=float),
                                                    np.arange(nkeep, dtype=float),
                                                    sources['xcentroid'][keep][srt],
                                                    sources['ycentroid'][keep][srt],
                                                    bkg[keep][srt],
                                                    rkron[keep][srt],
                                                    kflux_raw[keep][srt],
                                                    kflux[keep][srt],
                                                    np.ones(nkeep, dtype=float))),
                                axis=0)
        
        logger.info(f'Working on image {images[index]} ({len(stamps[keep])} stamps)')

        # Plot and save the stamps to a PDF file
        pdf_file = output_base.with_suffix(f".{images[index].path.stem.split('_')[0]}.stamps.pdf")
        show_stamps(all_stamps[index], pdf_file)

    # Save the source data to a text file
    header_base = 'Mean per-index Kron Radii:\n'
    for index in range(num_images):
        header_base += f'Image {images[index]}: {mean_rkron[index]:.2f}\n'
    np.savetxt(str(src_ofile), source_data,
                  fmt='%6.0f %4.0f %7.2f %7.2f %8.2f %7.2f %11.4e %11.4e %2.0f',
                  header=header_base + f"{'CHIP':>4} {'ID':>4} {'X':>7} {'Y':>7} {'BKG':>8} "
                                       f"{'R':>7} {'C_RAW':>11} {'C':>11} {'F':>2}")

    # Save the image data and stamps to a FITS file
    hdr = fits.Header()
    stamp_hdus = []
    for index in range(num_images):
        hdr[f'MRK_{index + 1:>02}'] = mean_rkron[index]
        if all_stamps[index] is None:
            continue
        stamp_hdus += [fits.ImageHDU(data=all_stamps[index].data, name=f'STAMPS_{index + 1:>02}')]
    fits.HDUList([fits.PrimaryHDU(header=hdr),
                  fits.ImageHDU(data=masked_images.data, name='IMAGES'),
                  fits.ImageHDU(data=np.ma.getmaskarray(masked_images).astype(np.uint8),
                                name='MASKS')]
                  + stamp_hdus).writeto(str(ofits), overwrite=True)
    
    return source_data

def main():

    return

if __name__ == '__main__':
    main()
