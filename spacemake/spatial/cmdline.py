from ..config import ConfigFile
from ..project_df import ProjectDF

from .he_integration import create_grayscale_expression_img
from ..util import message_aggregation, bool_in_str, str2bool

import argparse
import logging

logger_name = "spacemake.spatial"
logger = logging.getLogger(logger_name)

def get_expression_img_parser(with_umi_cutoff = False):
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)

    parser.add_argument('--project_id', type=str,
        required=True)

    parser.add_argument('--sample_id', type=str,
        required=True)

    parser.add_argument('--run_mode', type=str,
        required=True)

    if with_umi_cutoff:
        parser.add_argument('--umi_cutoff', type=int,
            required=True)

    parser.add_argument('--filter_percentage',
        type=int, required=False, default=70)

    return parser

def setup_spatial_parser(spmk, attach_to):
    parser = attach_to.add_parser('spatial',
        help = 'spacemake spatial commands')

    subparsers = parser.add_subparsers()

    grayscale_img_parser = subparsers.add_parser(
        'generate_grayscale_expression_img',
        parents=[get_expression_img_parser()])

    grayscale_img_parser.add_argument('--out_img',
        type=str,
        required=True)

    grayscale_img_parser.set_defaults(
        func=lambda args: generate_grayscale_expression_img_cmdline(spmk, args))

@message_aggregation(logger_name)
def generate_grayscale_expression_img_cmdline(spmk, args):
    import cv2 

    logger.info('Loading raw dge file...')
    adata_raw = spmk.load_raw_spatial_adata(
        project_id = args['project_id'],
        sample_id = args['sample_id'],
        run_mode_name = args['run_mode'])

    logger.info('Generating grayscale image...')
    img, img_be = create_grayscale_expression_img(adata_raw,
        filter_percentage=args['filter_percentage'])

    cv2.imwrite(args['out_img'], img)
