from ..config import ConfigFile
from ..project_df import ProjectDF

from .he_integration import create_aggregated_expression_img, \
    create_spot_expression_img
from ..util import message_aggregation, bool_in_str, str2bool
from ..errors import SpacemakeError

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

    parser.add_argument('--umi_cutoff', type=int,
        required=False)

    parser.add_argument('--binary_top_qth_percentile',
        type=int, required=False, default=30)

    parser.add_argument('--binary', type=str,
        required=False, default='False')

    parser.add_argument('--processed_data', type=str,
        required=False, default='False')

    parser.add_argument('--out_img',
        type=str,
        required=True)

    return parser

def setup_spatial_parser(spmk, attach_to):
    parser = attach_to.add_parser('spatial',
        help = 'spacemake spatial commands')

    subparsers = parser.add_subparsers()

    aggregated_img_parser = subparsers.add_parser(
        'create_aggregated_expression_img',
        parents=[get_expression_img_parser()])

    aggregated_img_parser.set_defaults(
        func=lambda args: create_expression_img_cmdline(spmk, args,
            'aggregated'))

    spot_img_parser = subparsers.add_parser(
        'create_spot_expression_img',
        parents=[get_expression_img_parser()])

    spot_img_parser.set_defaults(
        func=lambda args: create_expression_img_cmdline(spmk, args,
            'spot'))

@message_aggregation(logger_name)
def create_expression_img_cmdline(spmk, args, img_type):
    import cv2 

    logger.info('Loading dge file...')

    if str2bool(args['processed_data']):
        if not 'umi_cutoff' in args:
            raise SpacemakeError('When creating image from processed data,'
                ' a --umi_cutoff value must be provided') 

        adata = spmk.load_processed_adata(
            project_id = args['project_id'],
            sample_id = args['sample_id'],
            run_mode_name = args['run_mode'],
            umi_cutoff = args['umi_cutoff'])

    else:
        adata = spmk.load_raw_spatial_adata(
            project_id = args['project_id'],
            sample_id = args['sample_id'],
            run_mode_name = args['run_mode'])

    logger.info(f'Generating {img_type} expression image...')
    if img_type == 'spot':
        img, img_bw = create_spot_expression_img(adata,
            binary=str2bool(args['binary']))
    elif img_type == 'aggregated':
        img, img_bw = create_aggregated_expression_img(
            adata,
            binary_top_qth_percentile=int(args['binary_top_qth_percentile']))

    if str2bool(args['binary']):
        img = img_bw

    cv2.imwrite(args['out_img'], img)
