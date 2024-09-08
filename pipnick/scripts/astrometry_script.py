"""
Perform astrometric calibration on reduced images
""" 

from pipnick.scripts import scriptbase
from pipnick.utils.log import adjust_global_logger
import logging

class AstrometryPipeline(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Performs astrometric calibration on reduced images', width=width)
        parser.add_argument('maindir', type=str,
                            help='Path to the main directory containing the reduced directory with the FITS files to be astrometrically calibrated.')
        parser.add_argument('apikey', type=str,
                            help='API key from https://nova.astrometry.net account')
        parser.add_argument('-t', '--use_table', action='store_true',
                            help='Whether to use the table file to automatically exclude files that have been commented-out')
        parser.add_argument('-r', '--resolve', action='store_true', 
                            help="re-solves images with previously generated local solves")
        parser.add_argument('--excl_files', default=[], type=list,
                            help='List of file stems substrings to exclude (exact match not necessary).')
        parser.add_argument('--excl_objs', default=[], type=list,
                            help='List of object substrings to exclude (exact match not necessary).')
        parser.add_argument('--excl_filts', default=[], type=list,
                            help='List of filter substrings to exclude (exact match not necessary).')
        parser.add_argument('-vv', '--very_verbose', action='store_true', 
                            help="Display most detailed logs (use --verbosity for finer control)")
        parser.add_argument('--verbosity', default=4, type=str,
                            help='Level of verbosity to display (5=highest); overrides --verbose', 
                            choices=[1,2,3,4,5])
        return parser

    @staticmethod
    def main(args):
        
        from pipnick.pipelines.astrometry import astrometry_all
        
        if args.very_verbose:
            args.verbosity = 5
        log_levels = {1:'CRITICAL', 2:'ERROR', 3:'WARNING', 4:'INFO', 5:'DEBUG'}
        adjust_global_logger(log_levels[args.verbosity], __name__)
        logger = logging.getLogger(__name__)
        
        calib_files = astrometry_all(args.maindir, args.apikey, args.use_table, 
                                     args.resolve, args.excl_files, args.excl_objs, 
                                     args.excl_filts)
        
        return calib_files
        
        