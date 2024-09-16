import logging

logger_name = "spacemake.main"
logger = logging.getLogger(logger_name)


class Spacemake:
    """Spacemake.

    Class to access spacemake processed data from python.

    """

    def __init__(self, root):
        """__init__ constructor function of the Spacemake class

        :param root: Path to the spacemake root directory.
        :type root: str
        """
        from spacemake.config import get_global_config
        from spacemake.project_df import get_global_ProjectDF

        self.root = root
        self.config = get_global_config(root)
        self.project_df = get_global_ProjectDF(root)

    def load_processed_adata(
        self, project_id, sample_id, run_mode_name, umi_cutoff
    ): #-> anndata.AnnData:
        """Load spacemake processed data.

        :param project_id: project_id of the data to be loaded.
        :type project_id: str
        :param sample_id: sample_id of the data to be loaded.
        :type sample_id: str
        :param run_mode_name: name of the run mode of the data to be loaded.
            Each sample can have several run_modes during sample addition,
            here only one option needs to be provided.
        :type run_mode_name: str
        :param umi_cutoff: the umi_cutoff of the data to be loaded. Each
            run_mode can have several umi_cutoffs provided during configuration
            here only one option needs to be provided.
        :type umi_cutoff: int
        :returns: A spacemake processed and analyzed AnnData object, containing
            the results of the analysis.
        :rtype: anndata.AnnData
        """
        import scanpy as sc
        # import anndata


        self.project_df.assert_run_mode(project_id, sample_id, run_mode_name)
        run_mode = self.config.get_run_mode(run_mode_name)

        if not int(umi_cutoff) in [int(uc) for uc in run_mode.variables["umi_cutoff"]]:
            raise SpacemakeError(
                f"run_mode={run_mode} has no " + f"umi_cutoff={umi_cutoff}"
            )

        adata_raw = self.load_raw_spatial_adata(
            project_id=project_id, sample_id=sample_id, run_mode_name=run_mode_name
        )

        adata = sc.read(
            f"{self.root}/projects/{project_id}/processed_data/{sample_id}/"
            + f"illumina/complete_data/automated_analysis/{run_mode_name}/"
            + f"umi_cutoff_{umi_cutoff}/results.h5ad"
        )

        if "run_mode_variables" not in adata.uns.keys():
            adata.uns["run_mode_variables"] = run_mode.variables
        if "puck_variables" not in adata.uns.keys():
            adata.uns["puck_variables"] = adata_raw.uns["puck_variables"]

        return adata

    def load_raw_spatial_adata(
        self, project_id, sample_id, run_mode_name
    ): #-> anndata.AnnData:
        """Load raw, spacemake processed data.

        This function will load the raw countr matrix, created by spacemake.

        :param project_id: project_id of the raw data to be loaded.
        :type project_id: str
        :param sample_id: sample_id of the raw data to be loaded.
        :type sample_id: str
        :param run_mode_name: name of the run mode of the raw data to be loaded.
            Each sample can have several run_modes during sample addition,
            here only one option needs to be provided.
        :type run_mode_name: str
        :returns: A spacemake processed AnnData object, containing unfiltered
            raw expression data, and all cells or spatial units in the dataset.
        :rtype: anndata.AnnData
        """
        import scanpy as sc

        self.project_df.assert_run_mode(project_id, sample_id, run_mode_name)
        run_mode = self.config.get_run_mode(run_mode_name)

        dge_type = ""
        dge_cleaned = ""
        polyA_adapter_trimmed = ""
        mm_included = ""

        if run_mode.variables["polyA_adapter_trimming"]:
            polyA_adapter_trimmed = ".polyA_adapter_trimmed"

        if run_mode.variables["count_intronic_reads"]:
            dge_type = ".all"
        else:
            dge_type = ".exon"

        if run_mode.variables["count_mm_reads"]:
            mm_included = ".mm_included"

        if run_mode.variables["clean_dge"]:
            dge_cleaned = ".cleaned"

        adata = sc.read(
            f"{self.root}/projects/{project_id}/processed_data/{sample_id}/"
            + f"illumina/complete_data/dge/dge{dge_type}{dge_cleaned}"
            + f"{polyA_adapter_trimmed}{mm_included}.spatial_beads.h5ad"
        )

        if "puck_variables" not in adata.uns.keys():
            from spacemake.preprocess import attach_puck_variables

            adata = attach_puck_variables(
                adata,
                puck_variables=self.project_df.get_puck_variables(
                    project_id=project_id, sample_id=sample_id
                ),
            )

        if "run_mode_variables" not in adata.uns.keys():
            adata.uns["run_mode_variables"] = run_mode.variables

        return adata


def get_novosparc_variables(pdf, args):
    """get_novosparc_variables.

    :param pdf:
    :param args:
    """
    # assert that sample exists
    pdf.assert_sample(args["project_id"], args["sample_id"])

    def populate_variables_from_args(pdf, args, arg_prefix=""):
        """populate_variables_from_args.

        :param pdf:
        :param args:
        :param arg_prefix:
        """
        # get sample info
        sample_info = pdf.get_sample_info(
            project_id=args[f"{arg_prefix}project_id"],
            sample_id=args[f"{arg_prefix}sample_id"],
        )

        # populate return dictionary
        ret = {
            f"{arg_prefix}project_id": args[f"{arg_prefix}project_id"],
            f"{arg_prefix}sample_id": args[f"{arg_prefix}sample_id"],
        }

        # get run mode
        if f"{arg_prefix}run_mode" in args:
            ret[f"{arg_prefix}run_mode"] = args[f"{arg_prefix}run_mode"]
        else:
            run_mode_name = sample_info["run_mode"][0]
            ret[f"{arg_prefix}run_mode"] = run_mode_name
            logger.info(f"No run_mode provided, using {run_mode_name}")

        run_mode = pdf.config.get_run_mode(ret[f"{arg_prefix}run_mode"])

        if f"{arg_prefix}umi_cutoff" not in args:
            umi_cutoff = run_mode.variables["umi_cutoff"][0]
            ret[f"{arg_prefix}umi_cutoff"] = umi_cutoff
            logger.info(f"No umi_cutoff provided, using {umi_cutoff}")
        else:
            ret[f"{arg_prefix}umi_cutoff"] = args[f"{arg_prefix}umi_cutoff"]

        return ret

    ret = populate_variables_from_args(pdf, args)

    if "reference_project_id" not in args or "reference_sample_id" not in args:
        logger.info(
            "No reference_project_id or reference_sample_id provided,"
            + " running novosparc de-novo..."
        )
        ret["reference_project_id"] = ""
        ret["reference_sample_id"] = ""
        ret["reference_umi_cutoff"] = ""
        ret["reference_run_mode"] = ""
    else:
        pdf.assert_sample(args["reference_project_id"], args["reference_sample_id"])

        logger.info(
            "Using (project_id, sample_id)="
            + f"({args['reference_project_id']}, {args['reference_sample_id']})"
            + " reference, running novosparc with reference..."
        )

        novosparc_ret = populate_variables_from_args(pdf, args, arg_prefix="reference_")

        ret = {**ret, **novosparc_ret}

    return ret


_spacemake_instance = None


def get_spacemake_object():
    global _spacemake_instance
    if _spacemake_instance is None:
        _spacemake_instance = Spacemake(".")

    return _spacemake_instance


# def get_ConfigFile():
#     spmk = get_spacemake_object()
#     return spmk.config


# def get_ProjectDF():
#     spmk = get_spacemake_object()
#     return spmk.project_df
