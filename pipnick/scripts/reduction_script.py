"""
Perform reduction of raw astronomical data frames (overscan subtraction,
bias subtraction, flat division, cosmic ray masking)
""" 

from pipnick.scripts import scriptbase
from pipnick.convenience.log import adjust_global_logger
import logging

class ReductionPipeline(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Reduce images: subtract & trim overscan, subtract bias, divide flat', width=width)
        parser.add_argument('-dir', '--rawdir', default=None, type=str,
                            help='Directory with raw files to reduce.')
        parser.add_argument('-fin', '--table_path_in', default=None, type=str,
                            help='Path to input table file with raw FITS file information.')
        parser.add_argument('-fout', '--table_path_out', default='reduction_files_table.tbl', type=str,
                            help='Path to output table file for storing the raw FITS file information.')
        parser.add_argument('-s', '--save_inters', default=False, type=bool,
                            help='If True, save intermediate results during processing.')
        parser.add_argument('--excl_files', default=[], type=list,
                            help='List of file stems substrings to exclude (exact match not necessary).')
        parser.add_argument('--excl_obj_strs', default=[], type=list,
                            help='List of object substrings to exclude (exact match not necessary).')
        parser.add_argument('--excl_filts', default=[], type=list,
                            help='List of filter substrings to exclude (exact match not necessary).')
        parser.add_argument('-d', '--display', action='store_true', 
                            help="Display reduced images")
        parser.add_argument('-vv', '--very_verbose', action='store_true', 
                            help="Display most detailed logs (use --verbosity for finer control)")
        parser.add_argument('--verbosity', default=4, type=str,
                            help='Level of verbosity to display (5=highest); overrides --verbose', 
                            choices=[1,2,3,4,5])
        return parser

    @staticmethod
    def main(args):
        
        from pipnick.pipelines.reduction import reduce_all
        from pipnick.convenience.display_fits import display_many_nickel
        
        if args.very_verbose:
            args.verbosity = 5
        log_levels = {1:'CRITICAL', 2:'ERROR', 3:'WARNING', 4:'INFO', 5:'DEBUG'}
        adjust_global_logger(log_levels[args.verbosity], __name__)
        logger = logging.getLogger(__name__)
              
        logger.info("Running reduce_all()")
        red_files = reduce_all(args.rawdir, args.table_path_in, args.table_path_out,
                              args.save_inters, args.excl_files, args.excl_obj_strs, 
                              args.excl_filts)
        
        if args.display:
            display_many_nickel(red_files)
        
        