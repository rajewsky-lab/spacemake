import logging
import anndata

import numpy as np

from typing import List, Union
from spacemake.util import message_aggregation

logger_name = "spacemake.spatial.puck_collection"
logger = logging.getLogger(logger_name)
DEFAULT_REGEX_PUCK_ID = "(L[1-4][a-b]_tile_[1-2][0-7][0-9][0-9])"


def setup_parser(parser):
    """
    Set up command-line arguments for the script.

    :param parser: Argument parser object.
    :type parser: argparse.ArgumentParser
    :returns: Updated argument parser object.
    :rtype: argparse.ArgumentParser
    """
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
        default="outer",
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
    """
    Check if the observation keys in the AnnData object are unique.

    :param puck: AnnData object.
    :type puck: anndata.AnnData
    :param obs_key: Observation key to check uniqueness, defaults to "puck_id".
    :type obs_key: str, optional
    :returns: True if the observation keys are unique, False otherwise.
    :rtype: bool
    """
    return puck.obs[obs_key].nunique() == 1


def _transform_puck(puck: anndata.AnnData, pucks_transform: dict) -> anndata.AnnData:
    """
    Transform the spatial coordinates of the AnnData object.

    :param puck: AnnData object.
    :type puck: anndata.AnnData
    :param pucks_transform: Dictionary containing puck transformations.
    :type pucks_transform: dict
    :returns: Transformed AnnData object.
    :rtype: anndata.AnnData
    """
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
) -> anndata.AnnData:
    """
    Create a puck collection with transformed spatial coordinates.

    :param puck: AnnData object.
    :type puck: anndata.AnnData
    :param puck_transform: Dictionary containing puck transformations.
    :type puck_transform: dict
    :param reset_index: Flag to reset the observation index, defaults to True.
    :type reset_index: bool, optional
    :param transform: Flag to transform the spatial coordinates, defaults to True.
    :type transform: bool, optional
    :returns: Puck collection with transformed spatial coordinates.
    :rtype: anndata.AnnData
    """
    if reset_index:
        puck.obs_names = (
            puck.obs_names.astype(str) + ":" + puck.obs["puck_id"].astype(str)
        )

    if transform:
        puck = _transform_puck(puck, puck_transform)

    return puck


def parse_puck_id_from_path(f: str, puck_id_regex: str = DEFAULT_REGEX_PUCK_ID) -> str:
    """
    Parse the puck ID from the file path using the specified regex.

    :param f: File path.
    :type f: str
    :param puck_id_regex: Regular expression to find the puck ID, defaults to DEFAULT_REGEX_PUCK_ID.
    :type puck_id_regex: str, optional
    :returns: Puck ID extracted from the file path.
    :rtype: str
    """
    import os
    import re

    bname = os.path.basename(f)
    puck_id = re.findall(rf"{puck_id_regex}", bname)

    if len(puck_id) > 1:
        logger.warn(
            "Found more than one puck_id in the path. First one (index 0) will be used."
        )

    puck_id = puck_id[0]

    return puck_id


def read_pucks_to_list(
    f: Union[str, List[str]],
    puck_id: Union[int, List[int], None] = None,
    puck_id_regex: str = DEFAULT_REGEX_PUCK_ID,
    puck_id_key: str = "puck_id",
) -> List:
    """
    Read the pucks from file(s) and return a list of AnnData objects.

    :param f: File path(s).
    :type f: Union[str, List[str]]
    :param puck_id: Puck ID(s), defaults to None.
    :type puck_id: Union[int, List[int], None], optional
    :param puck_id_regex: Regular expression to find the puck ID in file names, defaults to DEFAULT_REGEX_PUCK_ID.
    :type puck_id_regex: str, optional
    :param puck_id_key: Name of the variable where puck IDs are stored, defaults to "puck_id".
    :type puck_id_key: str, optional
    :returns: List of AnnData objects representing the pucks.
    :rtype: List
    """

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


def parse_puck_coordinate_system_file(f: str) -> dict:
    """
    Parse the puck coordinate system file and return a dictionary.

    :param f: File path of the coordinate system file.
    :type f: str
    :returns: Dictionary representing the puck coordinate system.
    :rtype: dict
    """

    import pandas as pd

    cs = pd.read_csv(f, sep="[,|\t]", engine="python")

    cs = cs.set_index("puck_id")

    cs = cs.loc[~cs.index.duplicated(keep="first")]

    return cs.to_dict(orient="dict")


def merge_pucks_to_collection(
    pucks: List[str],
    puck_id: List[str],
    puck_coordinates: str,
    puck_id_regex: str = None,
    puck_id_key: str = "puck_id",
    no_reset_index: bool = False,
    no_transform: bool = False,
    merge_output: str = "same",
    join_output: str = "outer",
) -> anndata.AnnData:
    """
    Merge multiple pucks into a single puck collection.

    :param pucks: List of puck file paths.
    :type pucks: List[str]
    :param puck_id: List of puck IDs corresponding to the input files.
    :type puck_id: List[str]
    :param puck_coordinates: Name of the puck collection.
    :type puck_coordinates: str
    :param puck_id_regex: Regular expression to find the puck ID in file names, defaults to None.
    :type puck_id_regex: str, optional
    :param puck_id_key: Name of the variable where puck IDs are stored, defaults to "puck_id".
    :type puck_id_key: str, optional
    :param no_reset_index: Flag to not reset the observation index, defaults to False.
    :type no_reset_index: bool, optional
    :param no_transform: Flag to not transform the spatial coordinates, defaults to False.
    :type no_transform: bool, optional
    :param merge_output: How to merge pucks, can be "same", "unique", "first", or "only", defaults to "same".
    :type merge_output: str, optional
    :param join_output: How to join pucks, can be "inner" or "outer", defaults to "inner".
    :type join_output: str, optional
    :returns: Merged puck collection as an AnnData object.
    :rtype: anndata.AnnData
    """

    puck_transform = parse_puck_coordinate_system_file(puck_coordinates)

    pucks_list = read_pucks_to_list(pucks, puck_id, puck_id_regex, puck_id_key)

    puck_collection_list = []

    for puck in pucks_list:
        puck_collection_list += [
            create_puck_collection(
                puck,
                puck_transform,
                ~no_reset_index,
                ~no_transform,
            )
        ]

    puck_collection = anndata.concat(
        puck_collection_list, merge=merge_output, join=join_output
    )

    puck_collection.uns = {np.unique(puck.obs[puck_id_key])[0]: puck.uns for puck in puck_collection_list}

    return puck_collection

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
    puck_collection = merge_pucks_to_collection(
        pucks=args.pucks,
        puck_id=args.puck_id,
        puck_coordinates=args.puck_coordinates,
        puck_id_regex=args.puck_id_regex,
        puck_id_key=args.puck_id_key,
        no_reset_index=args.no_reset_index,
        no_transform=args.no_transform,
        merge_output=args.merge_output,
        join_output=args.join_output,
    )

    puck_collection.write_h5ad(args.output)


if __name__ == "__main__":
    cmdline()
