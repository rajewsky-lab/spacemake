import logging
import anndata

import numpy as np

from typing import List, Union
from spacemake.util import message_aggregation

logger_name = "spacemake.spatial.puck_collection"
logger = logging.getLogger(logger_name)
DEFAULT_REGEX_PUCK_ID = "(L[1-4][a-b]_tile_[1-2][0-7][0-9][0-9])"


def setup_parser(parser):
    parser.add_argument(
        "--pucks",
        type=str,
        nargs="+",
        help="path to spatial.h5ad AnnData file, one per puck",
        required=True,
    )

    parser.add_argument(
        "--puck-id",
        type=str,
        nargs="+",
        help="list of puck id for the input files, same order as pucks."
        + "Must be specified when the filenames do not contain a puck id that can be parsed with --puck-id-regex",
        default=None,
    )

    parser.add_argument(
        "--puck-coordinates",
        type=str,
        help="name of puck collection",
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        help="path to output.h5ad AnnData file",
        required=True,
    )

    parser.add_argument(
        "--puck-id-regex",
        type=str,
        help="regex to find puck id in file names",
        default=DEFAULT_REGEX_PUCK_ID,
    )

    parser.add_argument(
        "--puck-id-key",
        type=str,
        help="name of .obs variable where puck id are (will be) stored",
        default="puck_id",
    )

    parser.add_argument(
        "--merge-output",
        type=str,
        help='how to merge pucks, can be "same", "unique", "first", "only"',
        choices=["same", "unique", "first", "only"],
        default="same",
    )

    parser.add_argument(
        "--join-output",
        type=str,
        help='how to join pucks, can be "inner", "outer"',
        choices=["inner", "outer"],
        default="inner",
    )

    parser.add_argument(
        "--no-reset-index",
        default=False,
        action="store_true",
        help="do not reset the obs_name index of the AnnData object as  'obs_name:<puck_id_key>'; keep original 'obs_name'",
    )

    parser.add_argument(
        "--no-transform",
        default=False,
        action="store_true",
        help="do not transform the spatial coordinates of the AnnData object",
    )

    return parser


def check_obs_unique(puck: anndata.AnnData, obs_key: str = "puck_id") -> bool:
    return puck.obs[obs_key].nunique() == 1


def _transform_puck(puck: anndata.AnnData, pucks_transform: dict):
    if not check_obs_unique(puck, "puck_id"):
        raise ValueError(f"puck_id exists in AnnData object but are not unique")

    _puck_id = np.unique(puck.obs["puck_id"])[0]
    x_ofs = pucks_transform["x_offset"][_puck_id]
    y_ofs = pucks_transform["y_offset"][_puck_id]
    puck.obsm["spatial"][:, 0] += x_ofs
    puck.obsm["spatial"][:, 1] += y_ofs

    return puck


def create_puck_collection(
    puck: anndata.AnnData,
    puck_transform: dict,
    reset_index: bool = True,
    transform: bool = True,
):
    if reset_index:
        puck.obs_names = (
            puck.obs_names.astype(str) + ":" + puck.obs["puck_id"].astype(str)
        )

    if transform:
        puck = _transform_puck(puck, puck_transform)

    return puck


def parse_puck_id_from_path(f: str, puck_id_regex: str = DEFAULT_REGEX_PUCK_ID):
    import os
    import re

    bname = os.path.basename(f)
    puck_id = re.findall(rf"{puck_id_regex}", bname)

    if len(puck_id) > 1:
        logger.warn("Found more than one puck_id in the path. First one (index 0) will be used.")
    
    puck_id = puck_id[0]

    return puck_id


def read_pucks_to_list(
    f: Union[str, List[str]],
    puck_id: Union[int, List[int], None] = None,
    puck_id_regex: str = DEFAULT_REGEX_PUCK_ID,
    puck_id_key: str = "puck_id",
):
    import scanpy as sc

    if type(f) is str:
        f = [f]

    if puck_id is not None and type(puck_id) is str:
        puck_id = [puck_id]
    elif type(puck_id) is list and len(puck_id) != len(f):
        raise ValueError(
            f"Dimensions for f ({len(f)}) and puck_id ({len(puck_id)}) are not compatible"
        )

    pucks = []

    for i, f in enumerate(f):
        _f_obj = sc.read_h5ad(f)

        if "spatial" not in _f_obj.obsm.keys():
            raise ValueError(f"Could not find valid .obsm['spatial'] data in {f}")
        if puck_id_key not in _f_obj.obs.keys():
            if puck_id is None:
                _puck_id = parse_puck_id_from_path(f)
                if len(_puck_id) == 0:
                    raise ValueError(
                        f"Could not find a puck_id from the filename {f} with the regular expression {puck_id_regex}"
                    )
            else:
                _puck_id = puck_id[i]

            _f_obj.obs[puck_id_key] = _puck_id

        if puck_id_key != "puck_id":
            _f_obj.obs["puck_id"] = _puck_id

        if not check_obs_unique(_f_obj, "puck_id"):
            raise ValueError(
                f"puck_id exist in AnnData object but are not unique for the puck in file {f}"
            )

        pucks.append(_f_obj)

    return pucks


def parse_puck_coordinate_system_file(f: str):
    import pandas as pd

    cs = pd.read_csv(f, sep="[,|\t]", engine="python")

    cs = cs.set_index("puck_id")

    cs = cs.loc[~cs.index.duplicated(keep="first")]

    return cs.to_dict(orient="dict")


@message_aggregation(logger_name)
def cmdline():
    """cmdline."""
    import argparse

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="stitching/transforming pucks into a common global coordinate system",
    )
    parser = setup_parser(parser)

    args = parser.parse_args()

    puck_transform = parse_puck_coordinate_system_file(args.puck_coordinates)

    pucks_list = read_pucks_to_list(
        args.pucks, args.puck_id, args.puck_id_regex, args.puck_id_key
    )

    puck_collection_list = []

    for puck in pucks_list:
        puck_collection_list += [
            create_puck_collection(
                puck,
                puck_transform,
                ~args.no_reset_index,
                ~args.no_transform,
            )
        ]

    puck_collection = anndata.concat(
        puck_collection_list, merge=args.merge_output, join=args.join_output
    )

    puck_collection.write_h5ad(args.output)


if __name__ == "__main__":
    cmdline()
