Basic Data Reduction
====================

This tutorial guides you through the process of reducing raw astronomical
data using the `reduce_all()` function. The reduction process includes
overscan subtraction, bias subtraction, flat fielding and cosmic ray masking.

Overview
--------

The reduction pipeline is designed to process raw images and produce
reduced images that can be used for further analysis. The pipeline
handles tasks such as subtracting overscan, subtracting bias frames,
dividing by flat fields, and masking cosmic rays. The results are stored
in a directory structure stored in the input directory's parent.


Using the Command-Line Script
-----------------------------
For command-line usage, you can run the `pipnick_reduction`
script, which provides flexible options for performing reduction.


**Basic Use**

To execute the script, use the following command:

.. code::

  pipnick_reduction <maindir> [options]

Replace `<maindir>` with the directory containing the 'raw' directory.
All raw FITS files to be calibrated must be in /maindir/raw/.

Once you have run this script the first time, an ASCII Astropy table of
all files in /maindir/raw/ will be saved for reference /maindir/reduction_files.tbl.
You can even 'comment out' files with a `#` to be ignored in future runs.

.. code::

  pipnick_reduction <maindir> -t [options]

Toggle on use_table with `-t` to ignore all 'commented-out' files in the table.

**Using Optional Arguments**

The script also accepts the following command-line arguments:

- `-t` or `--use_table` (flag, optional)
  Whether to use the table file to automatically exclude files that have been commented-out.

- `-s` or `--save_inters` (bool, default=False):
  If `True`, save intermediate results during processing.

- `--excl_files` (list, optional):
  List of file stem substrings to exclude (exact match not necessary).

- `--excl_objs` (list, optional):
  List of object substrings to exclude (exact match not necessary).

- `--excl_filts` (list, optional):
  List of filter substrings to exclude (exact match not necessary).

- `-d` or `--display` (flag, optional):
  Display the reduced images.

- `-vv` or `--very_verbose` (flag, optional):
  Enable the most detailed logging. This option overrides `--verbosity`.

- `--verbosity` (int, default=4):
  Set the verbosity level for logging (1=CRITICAL, 5=DEBUG).


For example:

.. code::

  pipnick_reduction 'path/to/data/' --save_inters True --excl_files d1113 --excl_filts B --display

This command processes the raw files in the specified directory, saves
intermediate files, excludes certain files, and displays the reduced images.



Using the Reduction Function
----------------------------

To reduce your raw data using the `reduce_all()` function, follow these steps:

**Basic Use**

1. **Import the Required Functions**

   Begin by importing the `reduce_all()` function and any relevant utility
   functions from the `pipnick` package.

   .. code:: python

      from pipnick.pipelines.reduction import reduce_all
      from pipnick.convenience.display_fits import display_many_nickel

2. **Initialize Logging**

   Set up a logger to capture the output from the function. The verbosity
   level can be set to 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'.
   Logs will be displayed where the code is run and saved to a `.log` file
   at the 'DEBUG' level.

   .. code:: python

      import logging
      from pipnick.convenience.log import adjust_global_logger

      adjust_global_logger('INFO', __name__)
      logger = logging.getLogger(__name__)

3. **Specify the Raw Image Directory**

  ``maindir`` is the directory containing the 'raw' directory. All raw
  FITS files to be calibrated must be in /maindir/raw/.

  All results of the ``pipnick`` pipeline will be saved to this directory:
  reduction products will be saved to a folder called /reduced/ in ``maindir``.
  If ``save_inters`` is set to True, intermediate products will be saved
  to /processing/ in ``maindir``.

   .. code:: python

      maindir = 'path/to/maindir/'

4. **Run the Reduction Pipeline**

   The `reduce_all()` function reduces all files in /maindir/raw/, excluding
   files with ``'d1113'`` in the name or with ``'B'`` filter. It saves
   intermediate files (overscan subtracted, bias subtracted).

  This call also creates an ascii Astropy table of all files in /maindir/raw/
  for reference at /maindir/reduction_files.tbl, commenting out any files that
  were excluded.

   .. code:: python

      redfiles = reduce_all(maindir, save_inters=True, 
                            excl_files=['d1113'], excl_filts=['B'])

5. **Manual Exclusion of Files**

   The table can be edited with a `#` to comment out files (e.g., bad flats)
   that should be ignored in subsequent calls to `reduce_all()`.

   This call uses this table to determine exclusions. It will exclude the
   same files as in the first call, and adds in an exclusion for all files
   with ``'109'`` in the object name. These exclusions will be propagated
   to the table. This call also does not save intermediate files.

   .. code:: python

      redfiles = reduce_all(maindir, use_table=True, 
                            save_inters=False, excl_objs=['109'])

6. **Display the Reduced Files**

   After reduction, the reduced images can be displayed using the
   `display_many_nickel()` function.

   .. code:: python

      display_many_nickel(redfiles)


Viewing Results
---------------

Reduced images can be viewed using `display_many_nickel()` or in DS9. Note
that reduction may not correct certain "bad columns," which could be saturated
or otherwise problematic. These columns are masked according to definitions in
`pipnick.convenience.nickel_data`.
