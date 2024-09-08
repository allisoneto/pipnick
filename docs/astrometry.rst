Astrometric Solution
====================

This tutorial demonstrates how to use the `astrometry_all()` function
from the `pipnick` package to astrometrically calibrate images
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
able to solve the image. If a solution is not found in 60 seconds,
the function will recheck solution progress after all other images have
been calibrated. If a solution still has not been found, you may
manually check your submitted images in astrometry.net.

The astrometry and photometry pipelines are independent and can usually
run simultaneously. However, occasional conflicts may arise if both
pipelines attempt to access the same file concurrently.


Using the Command-Line Script
-----------------------------

For command-line usage, you can run the `pipnick_astrometry`
script, which provides flexible options for obtaining an astrometric solution.

**Basic Use**

To execute the script, use the following command:

.. code::

   pipnick_astrometry <maindir> <apikey> [options]

Replace `<maindir>` with the path to the directory containing the
raw & reduced directories. The reduced directory contains images to be
calibrated, and must be at /maindir/reduced/. Replace `<apikey>` with your
astrometry.net API key.

Once you have run this script the first time, an ASCII Astropy table of
all files in /maindir/reduced/ will be saved for reference at 
/maindir/astrometry_files.tbl. You may 'comment out' files to be
ignored in future runs with a `#`.

.. code::

  pipnick_astrometry <maindir> <api_key> -t [options]

Toggle on use_table with `-t` to ignore all 'commented-out' files in the table.

**Using Optional Arguments**

The script also accepts the following command-line arguments:

- `-t` or `--use_table` (flag, optional)
  Whether to use the table file to automatically exclude files that have been commented-out.

- `-r` or `--resolve` (flag, optional)
  If specified, re-solves images with previously generated local solves.

- `--excl_files` (list, optional):
  List of file stem substrings to exclude (exact match not necessary).

- `--excl_objs` (list, optional):
  List of object substrings to exclude (exact match not necessary).

- `--excl_filts` (list, optional):
  List of filter substrings to exclude (exact match not necessary).

- `-vv` or `--very_verbose` (flag, optional)
  Enables the most detailed logging. Overrides verbosity settings provided
  via `--verbosity`.

- `--verbosity` (int, optional)
  Sets the level of verbosity for logging. Acceptable values are 1 (CRITICAL),
  2 (ERROR), 3 (WARNING), 4 (INFO, default), and 5 (DEBUG). Overrides `--verbose`.

For example:

.. code::

   pipnick_astrometry '/path/to/maindir/ your_api_key' -t -r --excl_files d1113



Using Functions in Python File / Jupyter Notebook
-------------------------------------------------

To begin using the `astrometry_all()` function, follow these steps:

1. **Import the Function**

   First, import the `astrometry_all()` function from the
   `pipnick` package.

   .. code:: python

      from pipnick.pipelines.astrometry import astrometry_all

2. **Initialize Logging**

   Set up a logger to capture output from the functions. You can
   configure the verbosity level to 'DEBUG', 'INFO', 'WARNING',
   'ERROR', or 'CRITICAL'. Logs are displayed in the terminal or
   console and are always saved to a `.log` file at the 'DEBUG' level.

   .. code:: python

      import logging
      from pipnick.convenience.log import adjust_global_logger

      adjust_global_logger('INFO', __name__)
      logger = logging.getLogger(__name__)

3. **Specify the Image Directory**

   Define the main directory containing all data. The reduced
   directory contains images to be calibrated, and must be at
   /maindir/reduced/, and the results will be saved in a directory named
   `/data/astrometric/`.

   .. code:: python

      reddir = 'path/to/data/reduced/'

4. **Obtain an API Key**

   To use the astrometry.net service, you need an API key. Register an
   account at https://nova.astrometry.net/ and obtain your key from
   the "My Profile" section of the dashboard.

   .. code:: python

      api_key = "exampleapikey"

5. **Run the Astrometry Pipeline**

   Use the `astrometry_all()` function to process the images. This
   call This call saves astrometric solutions as WCS header to
   maindir/astrometric/, outputs the paths to these headers, and
   skips any images with pre-existing solutions to save time.

   As with the reduction pipeline, this call also creates an ascii
   Astropy table, but of all files in /maindir/reduced/ at
   /maindir/astrometry_files.tbl, commenting out any files that
   were excluded.

   .. code:: python

      calib_files = astrometry_all(reddir, api_key)

6. **Manual Exclusion of Files**

   As in the reduction pipeline, in the table created at
   /maindir/astrometry_files.tbl, you can comment out files
   (such as an unsolveable image) with a ``'#'`` to be ignored in
   a future call to ``astrometry_all``. Manual exclusions can be
   provided here as well, and will be propagated to the Astropy file table.

   In this second call, parameters can be changed to rely on this table
   for exclusion information, additionally exclude any image taken with the
   ``'B'`` filter, and to re-solve all images, regardless of whether
   previously solved.

   .. code:: python

      calib_files = astrometry_all(maindir, api_key, use_table=True,
                                   resolve=True, excl_filts=['B'])

Viewing Results
---------------

The header saved by this function contains information about the WCS
solution, which can be viewed in most text-viewing software.
