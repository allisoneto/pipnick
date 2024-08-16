import pkg_resources
import json
import logging.config


def load_logging_config():
    """
    Loads a logging configuration from a JSON file within the package.

    Returns
    -------
    dict
        The logging configuration loaded from the JSON file.
    """
    # Load the JSON configuration from a file within the 'nickelpipeline.convenience' package
    with pkg_resources.resource_stream('nickelpipeline.convenience', 'logging_config.json') as f:
        config = json.load(f)
    return config

def adjust_global_logger(log_level='INFO', name='all_others'):
    """
    Adjusts the global logging configuration to use a specific log level and output file.

    Parameters
    ----------
    log_level : str, optional
        The logging level to set for the console handler (default is 'INFO').
    name : str, optional
        The base name for the log file (default is 'all_others').

    Returns
    -------
    None
    """
    # Load the JSON configuration for logging
    output_file = f"log_{name.split('.')[-1]}.log"
    
    # Adjust the configuration for the file name and log level
    with pkg_resources.resource_stream('nickelpipeline.convenience', 'logging_config.json') as f:
        config = json.load(f)
        config['handlers']['file']['filename'] = output_file
        config['handlers']['console']['level'] = log_level

    # Configure logging with the loaded configuration
    logging.config.dictConfig(config)

def log_astropy_table(table):
    """
    Formats an Astropy Table into a string for logging purposes.

    Parameters
    ----------
    table : astropy.table.Table
        The Astropy Table to be formatted.

    Returns
    -------
    str
        The formatted string representation of the table.
    """
    # Convert each line of the table's formatted output into a string
    print_table = ""
    for line in table.pformat_all():
        print_table += f"{line}\n"
    return print_table
