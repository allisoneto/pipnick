"""
Perform reduction of raw astronomical data frames (overscan subtraction,
bias subtraction, flat division, cosmic ray masking)
""" 

from pipnick.scripts import scriptbase

class ReductionPipeline(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Perform basic image processing: (1) subtract & '
                                                'trim overscan, (2) subtract bias, and (3) '
                                                'divide flat', width=width)
        parser.add_argument('rawdir', type=str,
                            help='Path to parent directory with the raw data to be reduced.')
        parser.add_argument('-o', '--rdxdir', type=str,
                            help='Path for the reduced data.  If not provided, set to the parent '
                                 'directory of rawdir.')
        parser.add_argument('-t', '--table', default=None,
                            help='File with table with files to process')
        parser.add_argument('-s', '--save', default=False, action='store_true',
                            help='If True, save intermediate results during processing.')

        parser.add_argument('--excl_files', default=None, type=list,
                            help='List of file stems substrings to exclude (exact match not necessary).')
        parser.add_argument('--excl_objs', default=None, type=list,
                            help='List of object substrings to exclude (exact match not necessary).')
        parser.add_argument('--excl_filts', default=None, type=list,
                            help='List of filter substrings to exclude (exact match not necessary).')

        parser.add_argument('-d', '--display', action='store_true', help="Display reduced images")
        parser.add_argument('-vv', '--very_verbose', action='store_true', 
                            help="Display most detailed logs (use --verbosity for finer control)")
        parser.add_argument('--verbosity', default=4, type=str,
                            help='Level of verbosity to display (5=highest); overrides --verbose', 
                            choices=[1,2,3,4,5])
        return parser

    @staticmethod
    def main(args):
        
        from pipnick.utils.log import adjust_global_logger
        from pipnick.pipelines.reduction import reduce_all
        from pipnick.utils.display_fits import display_many_nickel

        # Initialize logger
        if args.very_verbose:
            args.verbosity = 5
        log_levels = {1:'CRITICAL', 2:'ERROR', 3:'WARNING', 4:'INFO', 5:'DEBUG'}
        adjust_global_logger(log_levels[args.verbosity], __name__)

        # Reduce the images
        metadata = reduce_all(args.rawdir, rdxdir=args.rdxdir, table=args.table, save=args.save,
                              excl_files=args.excl_files, excl_objs=args.excl_objs, 
                              excl_filts=args.excl_filts)

#        # Display reduced images
#        if args.display:
#            display_many_nickel(red_files)
