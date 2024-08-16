from nickelpipeline.pipelines.reduction import reduce_all
from nickelpipeline.pipelines.astrometry import astrometry_all
from nickelpipeline.pipelines.photometry import photometry_all
from nickelpipeline.pipelines.final_calib import final_calib_all
from nickelpipeline.convenience.graphs import plot_sources

import logging
from nickelpipeline.convenience.log import adjust_global_logger
adjust_global_logger('INFO', __name__)
logger = logging.getLogger(__name__)

# Paste your own API Key
api_key = 'exampleapikey'

# Define directory containing raw images
rawdir = 'data_example/raw'

# Basic reduction
red_files = reduce_all(rawdir=rawdir, save_inters=True)
reddir = red_files[0].parent.parent

# Astrometric calibration
astro_calib_files = astrometry_all(reddir, api_key)

# Photometric calibration
src_catalog_paths = photometry_all(reddir, group=True, plot_final=True,
                                   plot_inters=False)

# Final calibration (convert pixel coordinates -> RA/Dec)
photodir = src_catalog_paths[0].parent.parent
astrodir = astro_calib_files[0].parent
astrophot_data_tables = final_calib_all(photodir, astrodir)

# Display images & annotate sources
for object, src_table_dict in astrophot_data_tables.items():
    for file_key, src_table in src_table_dict.items():
        plot_sources(src_table, given_fwhm=8.0, flux_name='flux_psf', scale=1.5)