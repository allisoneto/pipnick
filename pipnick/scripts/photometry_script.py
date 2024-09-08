"""
Perform photometric calibration on reduced images
"""

from pipnick.scripts import scriptbase
from pipnick.utils.log import adjust_global_logger
import logging

class PhotometryPipeline(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Extracts sources from images & stores data in a table', width=width)
        parser.add_argument('maindir', type=str,
                            help='Path to main directory containing reduced directory with the files to be photometrically analyzed.')
        parser.add_argument('-t', '--use_table', action='store_true',
                            help='Whether to use the table file to automatically exclude files the have been commented-out')
        parser.add_argument('--excl_files', default=[], type=list,
                            help='List of file stems substrings to exclude (exact match not necessary).')
        parser.add_argument('--excl_objs', default=[], type=list,
                            help='List of object substrings to exclude (exact match not necessary).')
        parser.add_argument('--excl_filts', default=[], type=list,
                            help='List of filter substrings to exclude (exact match not necessary).')
        parser.add_argument('-th', '--thresh', default=8.0, type=float,
                            help='Threshold for source detection = background std * thresh.')
        parser.add_argument('-g', '--group', action='store_true', 
                            help='Consolidates groups of sources detected together into one source')
        parser.add_argument('-m', '--mode', default='all', type=str,
                            help='Mode to run photutils PSFPhotometry. `all` recommended.',
                            choices=['all', 'new'])
        parser.add_argument('-f', '--fittype', default='circ', type=str,
                            help='Which type of Moffat fit to use for PSF photometry',
                            choices=['circ', 'ellip'])
        parser.add_argument('-pf', '--plot_final', action='store_true', 
                            help="Displays images with sources & flux labelled")
        parser.add_argument('-pi', '--plot_inters', action='store_true', 
                            help="Displays images with initial sources & source groups for inspection")
        parser.add_argument('-vv', '--very_verbose', action='store_true', 
                            help="Display most detailed logs (use --verbosity for finer control)")
        parser.add_argument('--verbosity', default=4, type=str,
                            help='Level of verbosity to display (5=highest); overrides --verbose', 
                            choices=[1,2,3,4,5])
        return parser

    @staticmethod
    def main(args):
        
        from pipnick.pipelines.photometry import photometry_all

        if args.very_verbose:
            args.verbosity = 5
        log_levels = {1:'CRITICAL', 2:'ERROR', 3:'WARNING', 4:'INFO', 5:'DEBUG'}
        adjust_global_logger(log_levels[args.verbosity], __name__)
        logger = logging.getLogger(__name__)
        
        src_catalogs = photometry_all(args.maindir, use_table=args.use_table,
                                      excl_files=args.excl_files,
                                      excl_objs=args.excl_objs, 
                                      excl_filts=args.excl_filts, 
                                      thresh=args.thresh, group=args.group, 
                                      mode=args.mode, fittype=args.fittype,
                                      plot_final=args.plot_final, 
                                      plot_inters=args.plot_inters)
        
        return src_catalogs
 