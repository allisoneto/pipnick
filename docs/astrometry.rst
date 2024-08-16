Astrometric Solution
====================

This tutorial demonstrates how to use the `astrometry_all()` function
from the `nickelpipeline` package to astrometrically calibrate images
using astrometry.net.

Overview
--------

The `astrometry_all()` function performs astrometric calibration on
images that have undergone basic reduction. The function saves the
results as a FITS image (with a WCS file for coordinate conversion) and
a `.corr` table with detected source data and errors. It returns the
paths to one file type, based on the user's specifications.

Note that astrometry.net may take a significant amount of time to
process images. If too few stars are present, the service may not be
able to solve the image. The function is set to time out after
approximately 4 minutes.

The astrometry and photometry pipelines are independent and can usually
run simultaneously. However, occasional conflicts may arise if both
pipelines attempt to access the same file concurrently.


Using the Command-Line Script
-----------------------------

For command-line usage, you can run the `nickelpipeline_astrometry`
script, which provides flexible options for obtaining an astrometric solution.

**Basic Use**

To execute the script, use the following command:

.. code::

   nickelpipeline_astrometry <reddir> <apikey> [options]

Replace `<reddir>` with the path to your directory of reduced images and 
`<apikey>` with your astrometry.net API key. 

You can also specify optional arguments such as `--output_dir`, 
`--output_type`, `--resolve`, `--very_verbose`, and `--verbosity` 
to customize the script’s behavior.

**Using Optional Arguments**

- `-out` or `--output_dir` (str, optional)
  Path to the directory where the calibrated images will be saved. If not 
  specified, the default is `/astrometric/` in the same directory as `reddir`.

- `-t` or `--output_type` (str, optional)
  Defines whether to return paths to the calibrated image or a source table.
  Options are `image` (default) or `corr`.

- `-r` or `--resolve` (flag, optional)
  If specified, re-solves images with previously generated local solves.

- `-vv` or `--very_verbose` (flag, optional)
  Enables the most detailed logging. Overrides verbosity settings provided
  via `--verbosity`.

- `--verbosity` (int, optional)
  Sets the level of verbosity for logging. Acceptable values are 1 (CRITICAL),
  2 (ERROR), 3 (WARNING), 4 (INFO, default), and 5 (DEBUG). Overrides `--verbose`.

For example:

.. code::

   nickelpipeline_astrometry '/path/to/reduced/images your_api_key' -out /path/to/output -t corr



Using Functions in Python File / Jupyter Notebook
-------------------------------------------------

To begin using the `astrometry_all()` function, follow these steps:

1. **Import the Function**

   First, import the `astrometry_all()` function from the
   `nickelpipeline` package.

   .. code:: python

      from nickelpipeline.pipelines.astrometry import astrometry_all

2. **Initialize Logging**

   Set up a logger to capture output from the functions. You can
   configure the verbosity level to 'DEBUG', 'INFO', 'WARNING',
   'ERROR', or 'CRITICAL'. Logs are displayed in the terminal or
   console and are always saved to a `.log` file at the 'DEBUG' level.

   .. code:: python

      import logging
      from nickelpipeline.convenience.log import adjust_global_logger

      adjust_global_logger('INFO', __name__)
      logger = logging.getLogger(__name__)

3. **Specify the Image Directory**

   Define the directory containing the images to be calibrated. This
   can be the entire `/reduced/` directory or a specific object
   directory. The function will process files from all subdirectories.

   By default, the results will be saved in a directory named
   `/data/astrometric/` at the same level as `/data/reduced/`.

   .. code:: python

      reddir = 'path/to/data/reduced/'

4. **Obtain an API Key**

   To use the astrometry.net service, you need an API key. Register an
   account at https://nova.astrometry.net/ and obtain your key from
   the "My Profile" section of the dashboard.

   .. code:: python

      api_key = "exampleapikey"

5. **Run the Astrometry Pipeline**

   Use the `astrometry_all()` function to process the images. By
   default, the function saves the results to `data/astrometric/`,
   outputs the paths to the calibrated FITS images with WCS, and skips
   images with pre-existing solutions to save time.

   .. code:: python

      calib_files = astrometry_all(reddir, api_key)

6. **Customizing Parameters**

   You can customize the parameters to specify an output directory,
   change the file path output, and resolve all images, regardless of
   whether they have been previously solved.

   .. code:: python

      calib_files = astrometry_all(reddir, api_key, output_dir='path/to/output', 
                                   mode='corr', resolve=True)

Viewing Results
---------------

Astrometrically calibrated images can be viewed using DS9, which
automatically converts pixel coordinates to RA/Dec coordinates. The
header of the FITS image contains information about the WCS solution.
