Point-Source and Aperture Photometry
====================================

This tutorial demonstrates how to use the `photometry_all()` function
to perform photometric analysis on images that have undergone basic reduction.
Point-source photometry is performed using photutils.

Overview
--------

The `photometry_all()` function extracts and analyzes sources from
images and stores the results in `.csv` tables. These tables contain
data such as source positions and fluxes. Fluxes are calculated using PSF
photometry as well as aperture photometry. The function can also detect
multiple sources in one group, but groups can be consolidated.

The astrometry and photometry pipelines are independent and can usually
run simultaneously. However, conflicts may arise if both pipelines attempt
to access the same file concurrently.


Using the Command-Line Script
-----------------------------

For command-line usage, you can run the `pipnick_photometry`
script, which provides flexible options for executing photometric analysis.

**Basic Use**

To execute the script, use the following command:

.. code::

   pipnick_photometry <reddir> [options]

Replace `<reddir>` with the path to your directory of reduced images.

**Using Optional Arguments**

- `-out` or `--output_dir` (str, optional)
  Path to the directory where the results will be saved. If not 
  specified, the default is `/photometric/` in the same directory as `reddir`.

- `-t` or `--thresh` (float, optional)
  Threshold for source detection. The default is 8.0.

- `-g` or `--group` (flag, optional)
  Consolidates groups of sources detected together into one source.

- `-f` or `--fittype` (str, optional)
  Type of Moffat fit to use. Options are 'circ' (default) or 'ellip'.

- `-pf` or `--plot_final` (flag, optional)
  Displays images with sources and fluxes labeled.

- `-pi` or `--plot_inters` (flag, optional)
  Displays images with initial sources and source groups to determine
  by inspection which groups are valid.

- `-vv` or `--very_verbose` (flag, optional)
  Enables the most detailed logging. Overrides verbosity settings provided via `--verbosity`.

- `--verbosity` (int, optional)
  Sets the level of verbosity for logging. Acceptable values are 1 (CRITICAL),
  2 (ERROR), 3 (WARNING), 4 (INFO, default), and 5 (DEBUG). Overrides `--verbose`.

For example:

.. code::

   pipnick_photometry 'path/to/reduced/images' -out /path/to/output -t 10.0 -g -f ellip -pf -pi


Using the Photometry Function
-----------------------------

The same code can be run in a Python module with the same functionality.
To run the `photometry_all()` function, follow these steps:

**Basic Use**

1. **Import the Function**

   First, import the `photometry_all()` function from the
   `pipnick` package.

   .. code:: python

      from pipnick.pipelines.photometry import photometry_all

2. **Initialize Logging**

   Set up a logger to capture output from the function. You can
   configure the verbosity level to 'DEBUG', 'INFO', 'WARNING',
   'ERROR', or 'CRITICAL'. Logs are displayed in the terminal or
   console and are always saved to a `.log` file at the 'DEBUG' level.

   .. code:: python

      import logging
      from pipnick.convenience.log import adjust_global_logger

      adjust_global_logger('INFO', __name__)
      logger = logging.getLogger(__name__)

3. **Specify the Image Directory**

   Define the directory containing the images to be analyzed. This
   should be the entire `/reduced/` directory or a specific object
   directory. The function will process files from all subdirectories.

   By default, the results will be saved in a directory named
   `/data/photometric/` at the same level as `/data/reduced/`.

   .. code:: python

      reddir = 'path/to/data/reduced/'

4. **Run the Photometry Pipeline**

   Use the `photometry_all()` function to process the images. By
   default, the function saves the results in the default directory,
   applies a circular Moffat fit, and uses a detection threshold of
   8.0 times the background standard deviation.

   .. code:: python

      src_catalog_paths = photometry_all(reddir)

5. **Customizing Parameters**

   You can customize the function's behavior with various parameters.
   For example, you can set a different output directory, use an
   elliptical Moffat fit, consolidate source groups, or generate
   Matplotlib plots.

   .. code:: python

      src_catalog_paths = photometry_all(reddir, output_dir='path/to/output',
                                         thresh=15.0, group=True, fittype='ellip',
                                         plot_final=True, plot_inters=True)

Viewing Results
---------------

The output `.csv` files contain tables of detected sources with their
positions and fluxes. These tables are organized by object name and
saved in the specified output directory. If plotting options were
enabled, Matplotlib plots will show detected sources and source groups
for further inspection.