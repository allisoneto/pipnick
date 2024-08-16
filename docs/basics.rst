Basic Data Reduction
====================

Basic data-reduction steps include overscan subtraction, trimming,
bias-subtraction, and field-flattening.

To perform these steps, ...

Here are some options for how to run the code...

Image Reduction Pipeline
========================

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
For command-line usage, you can run the `nickelpipeline_reduction`
script, which provides flexible options for performing reduction.


**Basic Use**

To execute the script, use the following command:

.. code::

  nickelpipeline_reduction <rawdir> [options]

Replace `<rawdir>` with the path to your directory of raw images.

Once you have run this script the first time, an ASCII Astropy table of
all files in `rawdir` will be saved for reference. You can provide
the path to this table instead of the raw directory, and can even 'comment
out' files to be ignored.

.. code::

  nickelpipeline_reduction <table_path_in> [options]

Replace `<table_path_in>` with the path to this ASCII Astropy table.

**Using Optional Arguments**

The script accepts the following command-line arguments:

- `-dir` or `--rawdir` (str, optional):
  Directory containing raw files to reduce.

- `-fin` or `--table_path_in` (str, optional):
  Path to an input table file with raw FITS file information.

- `-fout` or `--table_path_out` (str, default='reduction_files_table.tbl'):
  Path to the output table file for storing raw FITS file information.

- `-s` or `--save_inters` (bool, default=False):
  If `True`, save intermediate results during processing.

- `--excl_files` (list, optional):
  List of file stem substrings to exclude (exact match not necessary).

- `--excl_obj_strs` (list, optional):
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

  nickelpipeline_reduction --rawdir 'path/to/data/raw/' --save_inters True --excl_files d1113 --excl_filts B --display

This command processes the raw files in the specified directory, saves intermediate files, excludes certain files, and displays the reduced images.



Using the Reduction Function
----------------------------

To reduce your raw data using the `reduce_all()` function, follow these steps:

**Basic Use**

1. **Import the Required Functions**

   Begin by importing the `reduce_all()` function and any relevant utility
   functions from the `nickelpipeline` package.

   .. code:: python

      from nickelpipeline.pipelines.reduction import reduce_all
      from nickelpipeline.convenience.display_fits import display_many_nickel

2. **Initialize Logging**

   Set up a logger to capture the output from the function. The verbosity
   level can be set to 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'.
   Logs will be displayed where the code is run and saved to a `.log` file
   at the 'DEBUG' level.

   .. code:: python

      import logging
      from nickelpipeline.convenience.log import adjust_global_logger

      adjust_global_logger('INFO', __name__)
      logger = logging.getLogger(__name__)

3. **Specify the Raw Image Directory**

   Define the directory containing the raw images. The processed results
   will be saved to a `/reduced/` folder within the parent directory of
   `rawdir`. If intermediate files are to be saved, they will be stored
   in a `/processing/` folder.

   .. code:: python

      rawdir = 'path/to/data/raw/'

4. **Run the Reduction Pipeline**

   The `reduce_all()` function will process all files in `rawdir`,
   excluding any files containing `'d1113'` in their name or those taken
   with a `'B'` filter. Intermediate products, such as overscan and
   bias-subtracted files, will be saved.

   Additionally, an ASCII Astropy table of all files in `rawdir` will be
   created for reference. Files that were excluded will be commented out
   in the table, which is saved by default as `files_table.tbl` in the
   parent directory of `rawdir`.

   .. code:: python

      redfiles = reduce_all(rawdir=rawdir, save_inters=True, 
                            excl_files=['d1113'], excl_filts=['B'])

5. **Manual Exclusion of Files**

   The created table can be edited to comment out files (e.g., bad flats)
   that should be ignored in subsequent calls to `reduce_all()`. The
   updated table must then be passed as `table_path_in`, instead of
   `rawdir`. Manual exclusions can also be provided, but they will only
   be recorded in the Astropy table if `table_path_out` is specified.

   .. code:: python

      redfiles = reduce_all(table_path_in='test_data/reduction_files_table2.tbl', 
                            table_path_out='test_data/reduction_files_table.tbl', 
                            save_inters=False, excl_obj_strs=['109'])

6. **Display the Reduced Files**

   After reduction, the reduced images can be displayed using the
   `display_many_nickel()` function.

   .. code:: python

      display_many_nickel(redfiles)


Viewing Results
---------------

Reduced images can be viewed using `display_many_nickel()` or in DS9. Note that reduction may not correct certain "bad columns," which could be saturated or otherwise problematic. These columns are masked according to definitions in `nickelpipeline.convenience.nickel_data`.
