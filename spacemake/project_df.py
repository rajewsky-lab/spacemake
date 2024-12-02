# import pandas as pd
import os
import math

# import argparse
# import datetime
import logging
import time

from spacemake.config import Puck
from spacemake.errors import *

from spacemake.config import ConfigFile
from spacemake.util import message_aggregation, assert_file, str_to_list
from spacemake.snakemake.variables import puck_barcode_files_summary
from typing import List, Dict

logger_name = "spacemake.project_df"
logger = logging.getLogger(logger_name)


class ProjectDF:
    """
    ProjectDF: class responsible for managing spacemake projects.

    :param file_path: path to the project_df.csv file, where we will save the project list.
    :param config: config file object
    :type config: ConfigFile
    :param df: A pandas dataframe, containing one row per sample
    :type df: pd.DataFrame
    """

    logger = logging.getLogger("spacemake.project_df.ProjectDF")

    # default values of the project dataframe columns
    project_df_default_values = {
        "puck_barcode_file_id": ["no_spatial_data"],
        "sample_sheet": None,
        "species": None,
        "demux_barcode_mismatch": 1,
        "demux_dir": None,
        "basecalls_dir": None,
        "R1": None,
        "R2": None,
        "reads": None,
        "longreads": None,
        "longread_signature": None,
        "investigator": "unknown",
        "sequencing_date": "unknown",
        "experiment": "unknown",
        "puck_barcode_file": None,
        "run_mode": ["default"],
        "barcode_flavor": "default",
        "is_merged": False,
        "merged_from": [],
        "puck": "default",
        "dge": None,
        "map_strategy": {  # default mapping strategy changes depending on if we have rRNA or not
            True: "rRNA:bowtie2->genome:STAR",  # map in parallel to rRNA and genome (default so far)
            False: "genome:STAR",  # w/o rRNA, map to genome directly
        },
        "adapter_flavor": "default",
    }

    project_df_dtypes = {
        "puck_barcode_file_id": "object",
        "sample_sheet": "str",
        "species": "str",
        "demux_barcode_mismatch": "int64",
        "demux_dir": "str",
        "basecalls_dir": "str",
        "R1": "object",
        "R2": "object",
        "reads": "object",
        "longreads": "str",
        "longread_signature": "str",
        "investigator": "str",
        "sequencing_date": "str",
        "experiment": "str",
        "puck_barcode_file": "object",
        "run_mode": "object",
        "barcode_flavor": "str",
        "is_merged": "bool",
        "merged_from": "object",
        "puck": "str",
        "dge": "str",
    }

    def __init__(self, file_path, config: ConfigFile = None):
        """__init__.

        :param file_path: path to pandas data frame (saved as .csv)
        :param config: ConfigFile
        :type config: ConfigFile
        """
        self.file_path = file_path
        self.config = config
        import pandas as pd

        if os.path.isfile(file_path):
            attempts = 0
            df_loaded_flag = False

            while not df_loaded_flag:
                try:
                    self.df = pd.read_csv(
                        file_path,
                        index_col=["project_id", "sample_id"],
                        na_values=["None", "none"],
                        dtype=self.project_df_dtypes,
                    )
                    df_loaded_flag = True
                except pd.errors.EmptyDataError as e:
                    if attempts < 5:
                        # wait 5 seconds before retrying
                        self.logger.warning(
                            f"The file '{self.file_path}' seems to be empty. Trying again ({attempts}/5) in 5 seconds"
                        )
                        time.sleep(5)
                        attempts = attempts + 1
                        continue
                    else:
                        raise e

            if self.df.empty:
                self.create_empty_df()
            else:
                # 'fix' the dataframe if there are inconsistencies
                self.fix()
        else:
            self.create_empty_df()

        self.assert_valid()

        self.logger = logging.getLogger(logger_name)

    def assert_valid(self):
        """assert_valid.

        this function iterates over projects/samples in the project_df, and asserts
        whether the specified variables are in accordance with the configuration file,
        and whether the specified files (R1, R2, puck_barcode_files) exist at the
        specified locations.
        """

        project_df_column_to_config_varname = {
            "run_mode": "run_modes",
            "species": "species",
            "barcode_flavor": "barcode_flavors",
            "puck": "pucks",
        }
        if not hasattr(self, "df"):
            raise SystemExit(
                ValueError(
                    "The 'project_df' does not exist in the ProjectDF object or is empty"
                )
            )
        elif self.df.empty:
            logger.warning("The 'project_df' in the ProjectDF object is empty")
            self.create_empty_df()

        _unique_samples = self._check_unique_samples()
        if _unique_samples is not None:
            raise SampleAlreadyExistsError(_unique_samples)

        for index, row in self.df.iterrows():
            # check that sample sheet file exists
            if (row["sample_sheet"] is not None) and (
                not os.path.exists(row["sample_sheet"])
            ):
                raise SystemExit(
                    FileNotFoundError(
                        f"At {index}, the 'sample_sheet' file does not exist"
                    )
                )

            # checking that variables are what they're supposed to be (according to config file)
            for pdf_col, config_var in project_df_column_to_config_varname.items():
                if type(row[pdf_col]) is list:
                    for _it_row in row[pdf_col]:
                        self.config.assert_variable(f"{config_var}", _it_row)
                elif type(row[pdf_col]) is str:
                    self.config.assert_variable(f"{config_var}", row[pdf_col])

            # check that puck_barcode_file(s), R1 and R2 files exist
            for n_col in ["puck_barcode_file", "R1", "R2"]:
                if n_col in ["R1", "R2"] and row["is_merged"]:
                    continue

                if type(row[n_col]) is list:
                    for _it_row in row[n_col]:
                        if not os.path.exists(_it_row):
                            raise SystemExit(
                                FileNotFoundError(
                                    f"At {index}, the {n_col} file does not exist: '{os.path.abspath(_it_row)}'"
                                )
                            )
                elif type(row[n_col]) is str and row[n_col] != "":
                    if not os.path.exists(row[n_col]):
                        raise SystemExit(
                            FileNotFoundError(
                                f"At {index}, the {n_col} file does not exist: '{os.path.abspath(_it_row)}'"
                            )
                        )

            # check that pucks are specified only if puck is spatial (or puck_collection is enabled)
            _valid_puck_coordinate = True
            if row["puck_barcode_file"] is None:
                _valid_puck_coordinate = False
            elif type(row["puck_barcode_file"]) is list:
                if len(row["puck_barcode_file"]) == 0:
                    _valid_puck_coordinate = False
            elif (
                type(row["puck_barcode_file"]) is str and row["puck_barcode_file"] == ""
            ):
                _valid_puck_coordinate = False

            # check that merging variables are properly specified
            # otherwise, this throws 'MissingInputException' because it cannot find the bam files
            if row["is_merged"]:
                if len(row["merged_from"]) < 2:
                    raise SystemExit(
                        SpacemakeError(
                            f"At {index}, there is <2 samples under 'merged_from'"
                        )
                    )
                for merged_i in row["merged_from"]:
                    if (
                        not isinstance(merged_i, tuple)
                        and not isinstance(merged_i, list)
                    ) or len(merged_i) != 2:
                        raise SystemExit(
                            SpacemakeError(
                                f"At {index}, wrong format for 'merged_from'.\n"
                                + "It must be something like "
                                + "[('project', 'sample_a'), ('project', 'sample_b')]"
                            )
                        )
                    try:
                        self.assert_sample(merged_i[0], merged_i[1])
                    except ProjectSampleNotFoundError as e:
                        self.logger.error(f"Merging error at {index}")
                        raise e

            if not _valid_puck_coordinate:
                _puck_vars = self.get_puck_variables(
                    project_id=index[0], sample_id=index[1]
                )
                if _puck_vars.get("coordinate_system", "") != "":
                    raise SystemExit(
                        SpacemakeError(
                            f"At {index}, the selected puck '{row['puck']}' "
                            + "contains a coordinate_system "
                            + "but no 'puck_barcode_files' are specified"
                        )
                    )

    def _check_unique_samples(self):
        for index, _ in self.df.iterrows():
            if len(self.df[self.df.index.isin([index])]) > 1:
                return index

        return None

    def create_empty_df(self):
        import pandas as pd

        index = pd.MultiIndex(
            names=["project_id", "sample_id"], levels=[[], []], codes=[[], []]
        )
        self.df = pd.DataFrame(self.project_df_default_values, index=index)

        self.df = self.df.astype(self.project_df_dtypes)
        self.dump()

    def compute_max_barcode_mismatch(self, indices: List[str]) -> int:
        """compute_max_barcode_mismatch.

        :param indices: List of illumina I7 index barcodes
        :type indices: List[str]
        :return: the maximum mismatch to be allowed for this set of index barcodes
        :rtype: int
        """
        num_samples = len(indices)

        if num_samples == 1:
            return 4
        else:
            max_mismatch = 3
            for i in range(num_samples - 1):
                for j in range(i + 1, num_samples):
                    hd = self.hamming_distance(indices[i], indices[j])
                    max_mismatch = min(max_mismatch, math.ceil(hd / 2) - 1)

        return max_mismatch

    def hamming_distance(self, string1: str, string2: str) -> int:
        """Cacluate hamming distance between two strings

        :param string1:
        :type string1: str
        :param string2:
        :type string2: str
        :rtype: int
        """
        return sum(c1 != c2 for c1, c2 in zip(string1, string2))

    def find_barcode_file(self, puck_barcode_file_id: str) -> str:
        """Tries to find path of a barcode file, using the puck_barcode_file_id.

        :param puck_barcode_file_id: puck_barcode_file_id of the puck we are looking for.
        :type puck_barcode_file_id: str
        :return: path of the puck file, containing barcodes, or None
        :rtype: str
        """
        # first find directory of puck file

        # return none or the path of the file
        def get_barcode_file(path):
            """get_barcode_file.

            :param path:
            """
            if os.path.isfile(path):
                return path

            return None

        def find_dir(name, path):
            """find_dir.

            :param name:
            :param path:
            """
            for root, dirs, files in os.walk(path):
                if name in dirs:
                    return os.path.join(root, name)

        puck_dir = find_dir(puck_barcode_file_id, self.config.puck_data["root"])
        path = None

        if puck_dir is not None:
            # puck dir exists, look for barcode file pattern
            path = os.path.join(puck_dir, self.config.puck_data["barcode_file"])

            return get_barcode_file(path)
        else:
            return self.project_df_default_values["puck_barcode_file"]

    def consolidate_pucks_merged_samples(self):
        for index, row in self.df.iterrows():
            project_id, sample_id = index
            puck_ids = row["puck_barcode_file_id"]

            if (not row["is_merged"]) or (
                not self.is_spatial(project_id, sample_id, puck_ids)
            ):
                continue

            if len(puck_ids) >= 1:
                puck_ids = puck_ids[0]
            elif len(puck_ids):
                puck_ids = self.project_df_default_values["puck_barcode_file_id"]

            merged_from = row["merged_from"]
            puck_id_file = set()

            for sample_tuple in merged_from:
                pid = self.df.loc[sample_tuple]["puck_barcode_file_id"]
                pbf = self.df.loc[sample_tuple]["puck_barcode_file"]
                _tuple = [(id, bf) for id, bf in zip(pid, pbf)]

                puck_id_file.update([tuple(t) for t in _tuple])

            pid, pbf = list(zip(*list(puck_id_file)))
            self.df.loc[index, "puck_barcode_file_id"] = list(pid)
            self.df.loc[index, "puck_barcode_file"] = list(pbf)

    def update_project_df_barcode_matches(self, prealigned=False):
        from spacemake.snakemake.variables import (
            puck_count_barcode_matches_summary,
            puck_count_prealigned_barcode_matches_summary,
        )
        import pandas as pd

        if prealigned:
            _bc_file = puck_count_prealigned_barcode_matches_summary
        else:
            _bc_file = puck_count_barcode_matches_summary

        for index, row in self.df.iterrows():
            project_id, sample_id = index

            puck_ids = row["puck_barcode_file_id"]
            if len(puck_ids) >= 1:
                puck_ids = puck_ids[0]
            elif len(puck_ids):
                puck_ids = self.project_df_default_values["puck_barcode_file_id"]

            if (row["is_merged"] and prealigned) or (
                not self.is_spatial(project_id, sample_id, puck_ids)
            ):
                continue

            # check if barcodes have been filtered
            _f_barcodes_df = _bc_file.format(project_id=project_id, sample_id=sample_id)

            # project_df is only updated if a prealignment barcode matching file is found
            if os.path.exists(_f_barcodes_df):
                barcodes_df = pd.read_csv(_f_barcodes_df)

                if "pass_threshold" not in barcodes_df.columns:
                    continue

                above_threshold_mask = barcodes_df["pass_threshold"] == 1

                _puck_barcode_files = barcodes_df[above_threshold_mask][
                    "puck_barcode_file"
                ].values.tolist()
                _puck_barcode_files_id = barcodes_df[above_threshold_mask][
                    "puck_barcode_file_id"
                ].values.tolist()

                self.df.at[index, "puck_barcode_file"] = _puck_barcode_files
                self.df.at[index, "puck_barcode_file_id"] = _puck_barcode_files_id

    def get_sample_info(self, project_id: str, sample_id: str) -> Dict:
        """get_sample_info.

        :param project_id:
        :type project_id: str
        :param sample_id:
        :type sample_id: str
        :return: A dictionary containing all the values of a given sample, from the ProjectDF.
        :rtype: Dict
        """
        # returns sample info from the projects df
        self.assert_sample(project_id, sample_id)
        out_dict = self.df.loc[(project_id, sample_id)].to_dict()

        return out_dict

    def is_external(self, project_id: str, sample_id: str) -> bool:
        """is_external.

        :param project_id:
        :type project_id: str
        :param sample_id:
        :type sample_id: str
        :rtype: bool
        """
        self.assert_sample(project_id, sample_id)

        data = self.df.loc[
            (project_id, sample_id),
            [
                "R1",
                "R2",
                "reads",
                "basecalls_dir",
                "sample_sheet",
                "longreads",
                "dge",
                "is_merged",
            ],
        ]

        self.logger.debug(
            f"project_id={project_id} sample_id={sample_id} "
            + f"R1,2={bool(data.R1 and data.R2)} "
            + f"basecall={(data.basecalls_dir and data.sample_sheet)} "
            + f"longreads={bool(data.longreads)} dge={bool(data.dge)} "
            + f"is_merged={bool(data.is_merged)}"
        )
        if (
            data.R2  # R1 is optional (bulk samples)
            or (data.basecalls_dir and data.sample_sheet)
            or (data.longreads)
            or (data.reads)
            and not data.dge
            or data.is_merged
        ):
            return False
        elif data.dge:
            return True
        else:
            raise SpacemakeError(
                f"Sample with id (project_id, sample_id)="
                + f"({project_id}, {sample_id}) is invalid."
            )

    def has_dge(self, project_id: str, sample_id: str) -> bool:
        """Returns True if a has dge. for Pacbio only samples returns False.

        :param project_id:
        :type project_id: str
        :param sample_id:
        :type sample_id: str
        :rtype: bool
        """
        self.assert_sample(project_id, sample_id)

        data = self.df.loc[
            (project_id, sample_id),
            [
                "R1",
                "R2",
                "reads",
                "basecalls_dir",
                "sample_sheet",
                "longreads",
                "dge",
                "is_merged",
            ],
        ]

        if (
            data.is_merged
            or data.R2
            or (data.sample_sheet and data.basecalls_dir)
            or data.dge
            or data.reads
        ):
            return True
        elif data.longreads:
            return False
        else:
            raise SpacemakeError(
                f"Sample with id (project_id, sample_id)="
                + f"({project_id}, {sample_id}) is invalid."
            )

    def is_spatial(
        self, project_id: str, sample_id: str, puck_barcode_file_id: str
    ) -> bool:
        """Returns true if a sample with index (project_id, sample_id) is spatial,
        meaning that it has spatial barcodes attached. Or, if the puck_barcode_file_id
        is 'puck_collection', meaning that the same is necessarily spatial (transformed
        from local to global coordinates)

        :param project_id:
        :type project_id: str
        :param sample_id:
        :type sample_id: str
        :rtype: bool
        """
        self.assert_sample(project_id, sample_id)
        puck_barcode_file = self.get_puck_barcode_file(
            project_id=project_id,
            sample_id=sample_id,
            puck_barcode_file_id=puck_barcode_file_id,
        )

        if puck_barcode_file is not None or puck_barcode_file_id == "puck_collection":
            return True
        else:
            return False

    def get_default_map_strategy_for_species(self, species):
        have_rRNA = "rRNA" in self.config.variables["species"][species]
        map_strategy = self.project_df_default_values["map_strategy"][have_rRNA]
        # print(
        #     f"getting default for '{species}' have_rRNA={have_rRNA} -> {map_strategy}"
        # )
        return map_strategy

    def fix(self):
        import numpy as np
        import pandas as pd

        modified = False
        # convert types
        self.df = self.df.where(pd.notnull(self.df), None)
        self.df = self.df.replace({np.nan: None})
        # replacing NaN with None

        # rename puck_id to puck_barcode_file_id, for backward
        # compatibility
        self.df.rename(
            columns={"puck_id": "puck_barcode_file_id"},
            inplace=True,
        )

        # convert list values stored as string
        self.df.run_mode = self.df.run_mode.apply(str_to_list)
        self.df.merged_from = self.df.merged_from.apply(str_to_list)

        # convert R1/R2 to list, if they are stored as string
        self.df.R1 = self.df.R1.apply(str_to_list)
        self.df.R2 = self.df.R2.apply(str_to_list)

        self.df.puck_barcode_file_id = self.df.puck_barcode_file_id.apply(str_to_list)
        self.df.puck_barcode_file = self.df.puck_barcode_file.apply(str_to_list)

        project_list = []
        # required if upgrading from pre-longread tree
        if not "reads" in self.df.columns:
            modified = True
            self.df["reads"] = None

        if not "longreads" in self.df.columns:
            modified = True
            self.df["longreads"] = None

        if not "longread_signature" in self.df.columns:
            modified = True
            self.df["longread_signature"] = None

        # required if upgrading from pre-bowtie2/map-strategy tree
        if not "map_strategy" in self.df.columns:
            modified = True
            self.df["map_strategy"] = self.df["species"].apply(
                self.get_default_map_strategy_for_species
            )

        # validate and correct map_strategies
        from spacemake.map_strategy import validate_mapstr

        corrected_map_strategies = []
        for row in self.df.itertuples():
            corrected_map_strategies.append(
                validate_mapstr(
                    row.map_strategy, config=self.config, species=row.species
                )
            )
        self.df["map_strategy"] = corrected_map_strategies

        if not "adapter_flavor" in self.df.columns:
            self.df["adapter_flavor"] = self.project_df_default_values["adapter_flavor"]
            # print("added adapter-flavor!")
            modified = True

        if modified:
            self.logger.warning(
                f".fix() reported changes! Saving migrated project_df.csv to '{self.file_path}'"
            )
            # self.logger.warning(self.df)
            # self.logger.warning(self.df.columns)
            # self.logger.warning(self.df["adapter_flavor"])
            self.dump()

        # per row updates
        # first create a series of a
        for ix, row in self.df.iterrows():
            s = pd.Series(self.project_df_default_values)

            # update puck barcode file info
            # for samples which have shared barcodes, and this barcode info is
            # stored in a puck, in the config file, before the id was set to
            # the name of the puck, and the puck_barcode_file was set to None.
            # Here we populate the puck_barcode_file into the path to the actual
            # file so that no errors are caused downstream.
            if (
                row["puck_barcode_file"] is None
                and row["puck_barcode_file_id"] is not None
            ):
                if len(row["puck_barcode_file_id"]) > 1:
                    raise SpacemakeError(
                        "When no barcode file provided, there "
                        + "only should be one id available"
                    )

                pbf_id = row["puck_barcode_file_id"][0]
                if pbf_id not in self.project_df_default_values["puck_barcode_file_id"]:
                    puck = self.config.get_puck(pbf_id)

                    row["puck_barcode_file"] = [puck.variables["barcodes"]]

            s.update(row)
            s.name = row.name
            project_list.append(s)

        self.df = pd.concat(project_list, axis=1).T
        self.df.is_merged = self.df.is_merged.astype(bool)
        self.df.index.names = ["project_id", "sample_id"]

    def get_puck_barcode_file(
        self, project_id: str, sample_id: str, puck_barcode_file_id: str
    ) -> str:
        if (
            puck_barcode_file_id
            in self.project_df_default_values["puck_barcode_file_id"]
        ):
            # if sample is not spatial, or we request the non-spatial puck
            return None
        else:
            ids = self.get_metadata(
                "puck_barcode_file_id", sample_id=sample_id, project_id=project_id
            )

            puck_barcode_files = self.get_metadata(
                "puck_barcode_file", sample_id=sample_id, project_id=project_id
            )

            # if no puck_barcode_file is provided, it means that barcode
            # file has to be fetched from the puck itself
            if puck_barcode_files is not None:
                for pid, pbf in zip(ids, puck_barcode_files):
                    if pid == puck_barcode_file_id:
                        return pbf

                return None

    def get_puck_barcode_ids_and_files(
        self,
        project_id: str,
        sample_id: str,
    ) -> str:
        puck_barcode_file_ids = self.get_metadata(
            "puck_barcode_file_id", project_id=project_id, sample_id=sample_id
        )

        puck_barcode_files = self.get_metadata(
            "puck_barcode_file", project_id=project_id, sample_id=sample_id
        )

        out_puck_barcode_files = []
        out_puck_barcode_file_ids = []

        # return only id-file pairs, for which file is not none
        if puck_barcode_files is not None:
            for pbf_id, pbf in zip(puck_barcode_file_ids, puck_barcode_files):
                out_puck_barcode_files.append(pbf)
                out_puck_barcode_file_ids.append(pbf_id)

        return out_puck_barcode_file_ids, out_puck_barcode_files

    def get_matching_puck_barcode_file_ids(
        self,
        project_id: str,
        sample_id: str,
    ):
        import pandas as pd

        summary_file = puck_barcode_files_summary.format(
            project_id=project_id, sample_id=sample_id
        )

        if not os.path.isfile(summary_file):
            # print(f"looking for summary file: '{summary_file}'")
            return self.project_df_default_values["puck_barcode_file_id"]

        df = pd.read_csv(summary_file)

        # Comment the following line that restricts the analysis of some tiles.
        # All tiles provided will now be processed -- it's up to the user to
        # provide a meaningful list of tiles.
        #
        # df = df.loc[(df.n_matching > 500) & (df.matching_ratio > 0.1)]

        pdf_ids = df.puck_barcode_file_id.to_list()

        pdf_ids.append(self.project_df_default_values["puck_barcode_file_id"][0])

        return pdf_ids

    def get_puck_barcode_file_metrics(
        self,
        project_id: str,
        sample_id: str,
        puck_barcode_file_id: str,
    ):
        import numpy as np

        summary_file = puck_barcode_files_summary.format(
            project_id=project_id, sample_id=sample_id
        )

        if not os.path.isfile(summary_file):
            return None
        import pandas as pd

        df = pd.read_csv(summary_file)

        df_puck = df.loc[df.puck_barcode_file_id == puck_barcode_file_id]

        # we calculate the stats for the puck_collection here
        # the max and min coordinates are not in global system
        if puck_barcode_file_id == "puck_collection":

            def multi_func(functions):
                def f(col):
                    return functions[col.name](col)

                return f

            df_pc = df._get_numeric_data().apply(
                multi_func(
                    {
                        "x_pos_min_px": np.min,
                        "x_pos_max_px": np.max,
                        "y_pos_min_px": np.min,
                        "y_pos_max_px": np.max,
                        "n_barcodes": np.mean,
                        "n_matching": np.mean,
                        "matching_ratio": np.mean,
                        "px_by_um": np.mean,
                    }
                )
            )

            return df_pc.to_dict()

        if df_puck.empty:
            return None
        else:
            return df_puck.iloc[0].to_dict()

    def get_puck(self, project_id: str, sample_id: str, return_empty=False) -> Puck:
        """get_puck.

        :param project_id: project_id of a sample
        :type project_id: str
        :param sample_id: sample_id of a sample
        :type sample_id: str
        :param return_empty:
        :return: A Puck object containing puck object
        :rtype: Puck
        """
        puck_name = self.get_metadata(
            "puck", project_id=project_id, sample_id=sample_id
        )

        return self.config.get_puck(puck_name, return_empty=return_empty)

    def get_puck_variables(
        self, project_id: str, sample_id: str, return_empty=False
    ) -> Dict:
        """get_puck_variables.

        :param project_id: project_id of a sample
        :type project_id: str
        :param sample_id: sample_id of a sample
        :type sample_id: str
        :param return_empty:
        :return: A dictionary containing the puck variables of a given sample
        :rtype: Dict
        """
        puck_name = self.get_metadata(
            "puck", project_id=project_id, sample_id=sample_id
        )

        return self.config.get_puck(puck_name, return_empty=return_empty).variables

    def get_metadata(self, field, project_id=None, sample_id=None, **kwargs):
        """get_metadata.

        :param field:
        :param project_id:
        :param sample_id:
        :param kwargs:
        """
        df = self.df
        if sample_id is not None:
            df = df.query("sample_id == @sample_id")

        if project_id is not None:
            df = df.query("project_id == @project_id")

        for key, value in kwargs.items():
            df = df.loc[df.loc[:, key] == value]

        dl = df.loc[:, field].to_list()
        if len(dl):
            return dl[0]
        else:
            return ""

    def dump(self):
        """dump."""
        self.logger.debug(f"writing project_df.csv to '{self.file_path}")
        self.df.to_csv(self.file_path)

    def add_sample_sheet(self, sample_sheet_path, basecalls_dir):
        """add_sample_sheet.

        :param sample_sheet_path:
        :param basecalls_dir:
        """
        import pandas as pd

        with open(sample_sheet_path) as sample_sheet:
            ix = 0
            investigator = None
            sequencing_date = None

            for line in sample_sheet:
                line = line.strip("\n")
                if "Investigator" in line:
                    investigator = line.split(",")[1]
                if "Date" in line:
                    sequencing_date = line.split(",")[1]
                if "[Data]" in line:
                    # the counter ix stops here
                    break
                else:
                    ix = ix + 1

        # read everything after [Data]
        df = pd.read_csv(sample_sheet_path, skiprows=ix + 1)
        # rename columns
        to_rename = {
            "Sample_ID": "sample_id",
            "Sample_Name": "puck_barcode_file_id",
            "Sample_Project": "project_id",
            "Description": "experiment",
            "index": "index",
        }
        df.rename(
            columns=to_rename,
            inplace=True,
        )
        # select only renamed columns
        df = df[to_rename.values()]
        df["species"] = df["experiment"].str.split("_").str[-1]
        df["investigator"] = investigator
        df["sequencing_date"] = sequencing_date

        # rename columns
        df["basecalls_dir"] = basecalls_dir
        df["demux_barcode_mismatch"] = self.compute_max_barcode_mismatch(df["index"])
        df["sample_sheet"] = sample_sheet_path
        df["demux_dir"] = (
            df["sample_sheet"].str.split("/").str[-1].str.split(".").str[0]
        )
        df["puck_barcode_file"] = df.puck_barcode_file_id.apply(self.find_barcode_file)
        df.set_index(["project_id", "sample_id"], inplace=True)

        for ix, row in df.iterrows():
            self.add_update_sample(project_id=ix[0], sample_id=ix[1], **row.to_dict())

    def assert_index_value(self, index_value, index_level):
        if not isinstance(index_value, list):
            index_value = [index_value]

        ixs = self.df.index.get_level_values(index_level)

        for ixv in index_value:
            if ixv not in ixs:
                raise ProjectSampleNotFoundError(index_level, ixv)

    def sample_exists(self, project_id=None, sample_id=None):
        """sample_exists.

        :param project_id:
        :param sample_id:
        """
        if project_id is None or sample_id is None:
            raise Exception(
                f"you need to provide a sample_id and project_id to check if sample exists"
            )
        else:
            ix = (project_id, sample_id)

            return ix in self.df.index

    def assert_sample(self, project_id, sample_id):
        if not self.sample_exists(project_id, sample_id):
            raise ProjectSampleNotFoundError(
                "(project_id, sample_id)", (project_id, sample_id)
            )

    def assert_run_mode(self, project_id, sample_id, run_mode_name):
        variables = self.get_sample_info(project_id, sample_id)

        if run_mode_name not in variables["run_mode"]:
            raise SpacemakeError(
                f"(project_id, sample_id)=({project_id},"
                + f"{sample_id}) has no run_mode={run_mode_name}\n"
                + f'run_mode has to be one of {variables["run_mode"]}'
            )

    def add_update_sample(
        self,
        action=None,
        project_id=None,
        sample_id=None,
        R1=None,
        R2=None,
        reads=None,
        dge=None,
        longreads=None,
        longread_signature=None,
        sample_sheet=None,
        basecalls_dir=None,
        is_merged=False,
        return_series=False,
        map_strategy=None,
        puck_barcode_file=None,
        puck_barcode_file_id=None,
        **kwargs,
    ):
        """add_update_sample.

        :param action:
        :param project_id:
        :param sample_id:
        :param R1:
        :param R2:
        :param reads:
        :param dge:
        :param longreads:
        :param longread_signature:
        :param sample_sheet:
        :param basecalls_dir:
        :param is_merged:
        :param return_series:
        :param map_strategy:
        :param kwargs:
        """
        import pandas as pd

        ix = (project_id, sample_id)
        sample_exists = self.sample_exists(*ix)

        if action is None and sample_exists:
            action = "update"
        elif action is None and not sample_exists:
            action = "add"
        elif action == "add" and sample_exists:
            raise SampleAlreadyExistsError(ix)
        elif action == "update" and not sample_exists:
            raise ProjectSampleNotFoundError("(project_id, sample_id)", ix)

        if action == "add":
            self.logger.info(f"Adding sample {ix}")
        elif action == "update":
            self.logger.info(f"Updating sample {ix}")
        else:
            raise ValueError(f"Unknown action {action}")

        # check variables
        # First we check if R1 and R2 is present.
        # If not, we check for longreads. If provided, we also need
        #   --longread-signature
        # If longreads not provided, we try with basecalls_dir and sample_sheet
        #   (only used by add_sample_sheet command)
        # If those area also not provided, we try to add a simple dge

        if R1 == ["None"]:
            R1 = None

        if R2 == ["None"]:
            R2 = None

        if reads == "None":
            reads = None

        if action == "add" and (R2 is None) and (reads is None) and not is_merged:
            self.logger.info("R2 not provided, trying longreads")

            if not longreads:
                self.logger.info(
                    "longreads not provided, trying basecalls_dir and sample_sheet"
                )
                if basecalls_dir is None or sample_sheet is None:
                    self.logger.info(
                        "basecalls_dir or sample_sheet not provided, trying dge"
                    )
                    if not dge:
                        raise SpacemakeError(
                            "Neither R1,R2, longreads, basecalls_dir & "
                            + "sample_sheet, nor dge were provided.\n"
                            + "Some reads/data has to be provided"
                        )
            else:
                if not longread_signature:
                    raise SpacemakeError(
                        "adding longreads requires to set --longread-signature as well (e.g. dropseq, chromium, noUMI, default, visium, slideseq_bc14,...)"
                    )

        # assert files first
        # if R1 is not None:
        assert_file(R1, default_value=None, extension=".fastq.gz")

        # if R2 is not None:
        assert_file(R2, default_value=None, extension=".fastq.gz")

        assert_file(reads, default_value=None, extension=[".txt", ".csv", ".tsv"])
        assert_file(longreads, default_value=None, extension="all")
        assert_file(
            dge,
            default_value=None,
            extension=[".h5", ".csv", ".h5ad", ".loom", ".txt", ".txt.gz"],
        )

        # assign reads
        # if there are strings, make them lists, as one sample can have many read files
        if R1 is not None and isinstance(R1, str):
            R1 = [R1]

        if R2 is not None and isinstance(R2, str):
            R2 = [R2]

        if R1 is not None and R2 is not None:
            if len(R1) != len(R2):
                raise SpacemakeError(
                    f"Trying to set an unmatching number of "
                    + f"read pairs for sample: {ix}.\n"
                    + f"# of R1 files = {len(R1)}\n"
                    + f"# of R2 files = {len(R2)}"
                )

        if "run_mode" in kwargs and isinstance(kwargs["run_mode"], str):
            # if a single run mode is provided as a string
            # create a list manually
            kwargs["run_mode"] = [kwargs["run_mode"]]

        # check if run mode exists
        for run_mode in kwargs.get("run_mode", []):
            if not self.config.variable_exists("run_modes", run_mode):
                raise ConfigVariableNotFoundError("run_modes", run_mode)

        config_variables_to_check = {
            "pucks": "puck",
            "barcode_flavors": "barcode_flavor",
            "adapter_flavors": "adapter_flavor",
            "species": "species",
        }

        for cv_plural, cv_singular in config_variables_to_check.items():
            if cv_singular in kwargs.keys():
                if not self.config.variable_exists(cv_plural, kwargs[cv_singular]):
                    raise ConfigVariableNotFoundError(cv_singular, kwargs[cv_singular])

        # if everything correct, add or update
        # first populate kwargs
        kwargs["R1"] = R1
        kwargs["R2"] = R2
        kwargs["reads"] = reads
        kwargs["dge"] = dge
        kwargs["longreads"] = longreads
        kwargs["longread_signature"] = longread_signature
        kwargs["sample_sheet"] = sample_sheet
        kwargs["basecalls_dir"] = basecalls_dir
        kwargs["is_merged"] = is_merged

        if action == "add" and map_strategy is None:
            # was not specified! Let's evaluate the default for this species
            # (depends on having rRNA reference or not)
            map_strategy = self.get_default_map_strategy_for_species(kwargs["species"])

        if (map_strategy is not None) and (map_strategy.endswith("-")):
            self.logger.warning(
                (
                    f"\n!!!!WARNING!!!\n map_strategy '{map_strategy}' ends with a trailing dash."
                    " Please ensure map_strategy is escaped with double-quotes, or the shell"
                    " interpretes the chevron in '->' as a redirect."
                )
            )

        if action == "add":
            _i_species = kwargs["species"]
        elif action == "update":
            _i_species = self.df.loc[ix]["species"]

        # validate and correct map_strategies
        from spacemake.map_strategy import validate_mapstr

        kwargs["map_strategy"] = validate_mapstr(
            map_strategy, config=self.config, species=_i_species
        )

        # TODO: remove
        # kwargs["map_strategy"] = map_strategy

        # populate puck_barcode_file
        if puck_barcode_file is not None:
            if isinstance(puck_barcode_file, str):
                puck_barcode_file = [puck_barcode_file]

            # if there are duplicates, raise error
            if len(puck_barcode_file) != len(set(puck_barcode_file)):
                raise SpacemakeError(
                    "Duplicate files provided for "
                    + "--puck_barcode_file. \n"
                    + f"files provided: {puck_barcode_file}"
                )

        if puck_barcode_file_id is not None:
            if isinstance(puck_barcode_file_id, str):
                puck_barcode_file_id = [puck_barcode_file_id]

            if len(puck_barcode_file_id) != len(set(puck_barcode_file_id)):
                raise SpacemakeError(
                    "Duplicate ids provided for "
                    + "--puck_barcode_file_id. \n"
                    + f"ids provided: {puck_barcode_file_id}"
                )

        # if puck barcode files and id's provided
        if puck_barcode_file_id is not None and puck_barcode_file is not None:
            # checklist
            # check if the lengths are the same
            if len(puck_barcode_file_id) != len(puck_barcode_file):
                raise SpacemakeError(
                    "Unmatching number of arguments provided"
                    + " for --puck_barcode_file and --puck_barcode_file_id.\n"
                    + "The provided number of elements should be the same"
                )

            kwargs["puck_barcode_file_id"] = puck_barcode_file_id
            kwargs["puck_barcode_file"] = puck_barcode_file
        elif puck_barcode_file_id is None and puck_barcode_file is not None:
            # if the user only provided puck_barcode_files
            self.logger.info(
                "No --puck_barcode_file_id provided. Generating" + " ids from filename"
            )

            puck_barcode_file_id = [
                os.path.basename(path).split(".")[0] for path in puck_barcode_file
            ]

            # if duplicates detected here: different files but same basename
            # this can happen if the file extensions are different
            if len(puck_barcode_file_id) != len(set(puck_barcode_file_id)):
                for i in range(len(puck_barcode_file_id)):
                    puck_barcode_file_id[i] = f"{puck_barcode_file_id[i]}_{i}"

            kwargs["puck_barcode_file_id"] = puck_barcode_file_id
            kwargs["puck_barcode_file"] = puck_barcode_file

        else:
            # if no puck barcode files are provided, we check if the puck has barcodes
            puck_name = kwargs.get("puck", None)

            if puck_name is not None:
                puck = self.config.get_puck(puck_name)

                if puck.has_barcodes:
                    kwargs["puck_barcode_file_id"] = [puck_name]
                    kwargs["puck_barcode_file"] = puck.variables["barcodes"]

        if sample_exists:
            new_project = self.df.loc[ix].copy()
            # pd.Series.update will only update values which are not None
            new_project.update(pd.Series(kwargs))
            self.df.loc[ix] = new_project
        else:
            new_project = pd.Series(self.project_df_default_values)
            new_project.name = ix
            new_project.update(kwargs)

            # after addition
            self.df = pd.concat([self.df, pd.DataFrame(new_project).T], axis=0)
            # as of pandas 2.0.1 the names of a MultiIndex do not survive concat.
            self.df.index.names = ["project_id", "sample_id"]

        if return_series:
            return (ix, new_project)

    def delete_sample(self, project_id, sample_id):
        """delete_sample.

        :param project_id:
        :param sample_id:
        """
        ix = (project_id, sample_id)

        if self.sample_exists(*ix):
            element = self.df.loc[ix]

            self.logger.info(
                f"Deleting sample: {ix}, with the following"
                + f" with the following variables:\n{element}"
            )

            self.df.drop(ix, inplace=True)
            return element
        else:
            raise ProjectSampleNotFoundError(
                "(project_id, sample_id)", (project_id, sample_id)
            )

    def add_samples_from_yaml(self, projects_yaml_file):
        """add_samples_from_yaml.

        :param projects_yaml_file:
        """
        import yaml

        config = yaml.load(open(projects_yaml_file), Loader=yaml.FullLoader)
        demux_projects = config.get("projects", None)

        if demux_projects is not None:
            # if we have projects in the config file
            # get the samples
            for ip in demux_projects:
                self.add_sample_sheet(ip["sample_sheet"], ip["basecalls_dir"])

        # add additional samples from config.yaml, which have already been demultiplexed.
        for project in config["additional_projects"]:
            self.add_update_sample(**project)

    def get_ix_from_project_sample_list(self, project_id_list=[], sample_id_list=[]):
        """get_ix_from_project_sample_list.

        :param project_id_list:
        :param sample_id_list:
        """

        # raise error if both lists are empty
        if project_id_list == [] and sample_id_list == []:
            raise NoProjectSampleProvidedError()

        # of only one provided use that, if both use intersection
        if project_id_list == []:
            ix = self.df.query("sample_id in @sample_id_list").index
        elif sample_id_list == []:
            ix = self.df.query("project_id in @project_id_list").index
        else:
            ix = self.df.query(
                "project_id in @project_id_list and sample_id in @sample_id_list"
            ).index

        return ix

    def set_variable(self, ix, variable_name, variable_key, keep_old=False):
        """set_variable.

        :param ix:
        :param variable_name:
        :param variable_key:
        :param keep_old:
        """
        variable_name_pl = self.config.main_variables_sg2pl[variable_name]
        self.config.assert_main_variable(variable_name_pl)
        self.config.assert_variable(variable_name_pl, variable_key)

        # get current value
        i_variable = self.df.at[ix, variable_name]

        # if list, we either append or not
        if isinstance(i_variable, list):
            if keep_old:
                # keep the old list as well
                i_variable = list(set(i_variable + variable_key))
            else:
                i_variable = list(set(variable_key))

        # if we do not keep the old, simply set
        else:
            i_variable = variable_key

        if i_variable == [] or i_variable is None:
            raise EmptyConfigVariableError(variable_name)

        self.df.at[ix, variable_name] = i_variable

    def remove_variable(self, ix, variable_name, variable_key):
        """remove_variable.

        :param ix:
        :param variable_name:
        :param variable_key:
        """
        variable_name_pl = self.config.main_variables_sg2pl[variable_name]
        self.config.assert_main_variable(variable_name_pl)
        self.config.assert_variable(variable_name_pl, variable_key)

        if not isinstance(variable_key, list):
            raise ValueError("variable_key has to be a list")

        i_variable = self.df.at[ix, variable_name]

        i_variable = [val for val in i_variable if val not in variable_key]
        if i_variable == [] or i_variable is None:
            raise EmptyConfigVariableError(variable_name)

        self.df.at[ix, variable_name] = i_variable

    def set_remove_variable(
        self,
        variable_name,
        variable_key,
        action,
        project_id_list=[],
        sample_id_list=[],
        keep_old=False,
    ):
        """set_remove_variable.

        :param variable_name:
        :param variable_key:
        :param action:
        :param project_id_list:
        :param sample_id_list:
        :param keep_old:
        """
        self.assert_projects_samples_exist(project_id_list, sample_id_list)

        variable_name_pl = self.config.main_variables_sg2pl[variable_name]
        self.config.assert_main_variable(variable_name_pl)

        ix = self.get_ix_from_project_sample_list(
            project_id_list=project_id_list, sample_id_list=sample_id_list
        )

        if isinstance(variable_key, list):
            for var in variable_key:
                self.config.assert_variable(variable_name_pl, var)
        else:
            self.config.assert_variable(variable_name_pl, variable_key)

        for i, row in self.df.loc[ix, :].iterrows():
            # add/remove variable. if it already exists dont do anything
            if action == "set":
                self.set_variable(i, variable_name, variable_key, keep_old)
            elif action == "remove":
                self.remove_variable(i, variable_name, variable_key)

        return ix.to_list()

    def assert_projects_samples_exist(self, project_id_list=[], sample_id_list=[]):
        """assert_projects_samples_exist.

        :param project_id_list:
        :param sample_id_list:
        """
        for project in project_id_list:
            if project not in self.df.index.get_level_values("project_id"):
                raise ProjectSampleNotFoundError("project_id", project)

        for sample in sample_id_list:
            if sample not in self.df.index.get_level_values("sample_id"):
                raise ProjectSampleNotFoundError("sample_id", sample)

    def merge_samples(
        self,
        merged_project_id,
        merged_sample_id,
        project_id_list=[],
        sample_id_list=[],
        **kwargs,
    ):
        """merge samples.

        :param merged_project_id:
        :param merged_sample_id:
        :param project_id_list:
        :param sample_id_list:
        :param kwargs:
        """
        # check if projects and samples with these IDs exist
        self.assert_projects_samples_exist(project_id_list, sample_id_list)

        ix = self.get_ix_from_project_sample_list(
            project_id_list=project_id_list, sample_id_list=sample_id_list
        )

        consistent_variables = list(self.config.main_variables_sg2pl.keys())
        consistent_variables.remove("run_mode")
        consistent_variables.remove("adapter")
        consistent_variables.remove("quant")

        # append variable "map-strategy" manually.
        # TODO: consider moving "map_strategy" into the main_variables_sg2pl, needs additional parser etc.
        consistent_variables.append("map_strategy")

        ix_list = ix.to_list()

        if ix_list == []:
            raise ProjectSampleNotFoundError(
                "(project_id_list, sample_id_list)", (project_id_list, sample_id_list)
            )

        # check for variable inconsistency
        # raise error if variable different between samples
        for variable in consistent_variables:
            if variable in kwargs and variable not in ["species", "barcode_flavor"]:
                # if variable provided from command line, skip
                self.logger.info(f"{variable} provided, skipping deduction...")
                continue

            self.logger.info(
                f"{variable} not provided. deducing from merged samples..."
            )

            variable_val = self.df.loc[ix, variable].to_list()

            if len(set(variable_val)) > 1:
                raise InconsistentVariablesDuringMerge(
                    variable_name=variable, variable_value=variable_val, ix=ix
                )
            else:
                # attach the deduced, consisten variable
                kwargs[variable] = variable_val[0]

        # get puck_barcode_files
        if "puck_barcode_file" not in kwargs or "puck_barcode_file_id" not in kwargs:
            kwargs["puck_barcode_file_id"] = []
            kwargs["puck_barcode_file"] = []

            for _, row in self.df.loc[ix].iterrows():
                if row["puck_barcode_file"] is None:
                    continue
                else:
                    for pbf_id, pbf in zip(
                        row["puck_barcode_file_id"], row["puck_barcode_file"]
                    ):
                        if (
                            pbf_id not in kwargs["puck_barcode_file_id"]
                            and pbf not in kwargs["puck_barcode_file"]
                        ):
                            kwargs["puck_barcode_file_id"].append(pbf_id)
                            kwargs["puck_barcode_file"].append(pbf)

        # after all checks, log that we are merging
        self.logger.info(f"Merging samples {ix_list} together\n")
        variables_to_deduce = ["investigator", "experiment", "sequencing_date"]

        for variable in variables_to_deduce:
            if variable not in kwargs:
                kwargs[variable] = ";".join(self.df.loc[ix, variable].unique())

        # if no run_mode provided, overwrite with user defined one
        if "run_mode" not in kwargs.keys():
            run_mode_lists = [
                set(run_mode) for run_mode in self.df.loc[ix].run_mode.to_list()
            ]

            # join run modes from parent samples
            if len(run_mode_lists) == 1:
                run_mode = run_mode_lists[0]
            else:
                run_mode = run_mode_lists[0].intersection(*run_mode_lists[1:])

            # create a list from the set intersection
            run_mode = list(run_mode)

            # if there are no common elements, throw an error
            if len(run_mode) == 0:
                raise InconsistentVariablesDuringMerge(
                    variable_name="run_mode", variable_value=run_mode_lists, ix=ix_list
                )

            # finally add run mode to arguments
            kwargs["run_mode"] = run_mode

        # set the action to add
        kwargs["action"] = "add"

        sample_added = self.add_update_sample(
            project_id=merged_project_id,
            sample_id=merged_sample_id,
            is_merged=True,
            merged_from=ix.to_list(),
            **kwargs,
        )

        return (sample_added, ix)


__global_ProjectDF = None


def get_global_ProjectDF(root="."):
    global __global_ProjectDF
    if __global_ProjectDF is None:
        from spacemake.config import get_global_config

        __global_ProjectDF = ProjectDF(
            f"{root}/project_df.csv", config=get_global_config()
        )

    return __global_ProjectDF
