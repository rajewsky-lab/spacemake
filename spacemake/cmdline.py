import argparse
import os
import logging
from spacemake.util import message_aggregation
from spacemake.config import ConfigFile
from spacemake.contrib import __version__, __license__, __author__, __email__

logger_name = "spacemake.main"
logger = logging.getLogger(logger_name)

def get_project_sample_parser(allow_multiple=False, prepend="", help_extra=""):
    """
    Return a parser for project_id's and sample_id's

    :param allow_multiple: if true, we allow multiple projects and samples,
        and the parser will have a `--project_id_list` and a `--sample_id_list`
        parameter
    :param prepend: a string that will be prepended before the arguments
    :param help_extra: extra help message
    :return: a parser object
    :rtype: argparse.ArgumentParser
    """
    logger.debug(f"get_project_sample_parser(prepend={prepend}) called")

    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)
    project_argument = "project_id"
    sample_argument = "sample_id"
    required = True
    nargs = None
    default = None

    if allow_multiple:
        project_argument = f"{project_argument}_list"
        sample_argument = f"{sample_argument}_list"
        required = False
        nargs = "*"
        default = []

    parser.add_argument(
        f"--{prepend}{project_argument.replace('_', '-')}",
        type=str,
        required=required,
        nargs=nargs,
        default=default,
        help=f"{project_argument} {help_extra}",
        dest=f"{prepend.replace('-', '_')}{project_argument}",
    )
    parser.add_argument(
        f"--{prepend}{sample_argument.replace('_', '-')}",
        type=str,
        required=required,
        nargs=nargs,
        default=default,
        help=f"{sample_argument} {help_extra}",
        dest=f"{prepend.replace('-', '_')}{sample_argument}",
    )

    # Add legacy snake_case arguments
    parser.add_argument(
        f"--{prepend}{project_argument}",
        type=str,
        required=False,  # Not required since kebab-case takes precedence
        nargs=nargs,
        default=None,  # No default here; rely on kebab-case
        dest=f"{prepend.replace('-', '_')}{project_argument}",
    )
    parser.add_argument(
        f"--{prepend}{sample_argument}",
        type=str,
        required=False,
        nargs=nargs,
        default=None,
        dest=f"{prepend.replace('-', '_')}{sample_argument}",
    )

    return parser


def get_add_sample_sheet_parser():
    """
    Returns parser for sample sheet addition

    :return: parser
    :rtype: argparse.ArgumentParser
    """
    logger.debug("get_add_sample_sheet_parser() called")
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="add a new sample sheet to the samples",
        add_help=False,
    )

    parser.add_argument(
        "--sample-sheet",
        type=str,
        help="the path to the Illumina sample sheet",
        required=True,
        dest="sample_sheet",
    )
    parser.add_argument(
        "--sample_sheet",
        type=str,
        help=argparse.SUPPRESS,
        required=False,
        dest="sample_sheet",
    )
    
    parser.add_argument(
        "--basecalls-dir",
        type=str,
        help="path to the basecalls directory",
        required=True,
        dest="basecalls_dir",
    )
    parser.add_argument(
        "--basecalls_dir",
        type=str,
        help=argparse.SUPPRESS,
        required=False,
        dest="basecalls_dir",
    )

    return parser


def get_sample_main_variables_parser(
    defaults=False,
    species_required=False,
    main_variables=[
        "barcode_flavor",
        "adapter_flavor",
        "species",
        "puck",
        "run_mode",
        "map_strategy",
    ],
):
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)
    logger.debug(
        f"get_sample_main_variables_parser called with main_variables={main_variables}"
    )

    if "barcode_flavor" in main_variables:
        parser.add_argument(
            "--barcode-flavor",
            type=str,
            default="default" if defaults else None,
            help="barcode flavor for this sample",
            dest="barcode_flavor",
        )
        parser.add_argument(
            "--barcode_flavor",
            type=str,
            default="default" if defaults else None,
            help=argparse.SUPPRESS,
            dest="barcode_flavor",
        )

    if "adapter_flavor" in main_variables:
        parser.add_argument(
            "--adapter-flavor",
            type=str,
            default="default" if defaults else None,
            help="barcode flavor for this sample",
            dest="adapter_flavor",
        )
        parser.add_argument(
            "--adapter_flavor",
            type=str,
            default="default" if defaults else None,
            help=argparse.SUPPRESS,
            dest="adapter_flavor",
        )

    if "species" in main_variables:
        parser.add_argument(
            "--species",
            type=str,
            help="add the name of the species for this sample",
            required=species_required,
        )

    if "map_strategy" in main_variables:
        parser.add_argument(
            "--map-strategy",
            type=str,
            help="string constant definining mapping strategy. Can be multi-stage and use bowtie2 or STAR (see documentation)",
            required=False,
            default="STAR:genome",
            dest="map_strategy",
        )
        parser.add_argument(
            "--map_strategy",
            type=str,
            help=argparse.SUPPRESS,
            required=False,
            default=None,
            dest="map_strategy",
        )

    if "puck" in main_variables:
        parser.add_argument(
            "--puck",
            type=str,
            help="name of the puck for this sample. if puck_type contains a \n"
            + "`barcodes` path to a coordinate file, those coordinates\n"
            + " will be used when processing this sample. if \n"
            + " not provided, a default puck will be used with \n"
            + "width_um=3000, spot_diameter_um=10",
        )

        parser.add_argument(
            "--puck-barcode-file-id",
            type=str,
            help="puck_barcode_file_id of the sample to be added/update",
            nargs="+",
            dest="puck_barcode_file_id",
        )
        parser.add_argument(
            "--puck_barcode_file_id",
            type=str,
            help=argparse.SUPPRESS,
            nargs="+",
            dest="puck_barcode_file_id",
        )

    parser.add_argument(
        "--puck-barcode-file",
        type=str,
        nargs="+",
        help="the path to the file contining (x,y) positions of the barcodes",
        dest="puck_barcode_file",
    )
    parser.add_argument(
        "--puck_barcode_file",
        type=str,
        nargs="+",
        help=argparse.SUPPRESS,
        dest="puck_barcode_file",
    )

    if "run_mode" in main_variables:
        parser.add_argument(
            "--run-mode",
            type=str,
            nargs="+",
            help="run_mode names for this sample.\n"
            + "the sample will be processed using the provided run_modes.\n"
            + "for merged samples, if left empty, the run_modes of the \n"
            + "merged (input) samples will be intersected.\n",
            dest="run_mode",
        )
        parser.add_argument(
            "--run_mode",
            type=str,
            nargs="+",
            help=argparse.SUPPRESS,
            dest="run_mode",
        )

    return parser


def get_sample_extra_info_parser():
    logger.debug("get_sample_extra_info_parser() called")
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)

    parser.add_argument(
        "--investigator",
        type=str,
        help="add the name of the investigator(s) responsible for this sample",
    )

    parser.add_argument("--experiment", type=str, help="description of the experiment")

    import datetime

    parser.add_argument(
        "--sequencing-date",
        type=datetime.date.fromisoformat,
        help="sequencing date of the sample. format: YYYY-MM-DD",
        dest="sequencing_date",
    )
    parser.add_argument(
        "--sequencing_date",
        type=datetime.date.fromisoformat,
        help=argparse.SUPPRESS,
        dest="sequencing_date",
    )

    return parser


def get_data_parser(reads_required=False):
    """
    Returns a parser which contain extra arguments for a given sample.
    The returned parser will contain the --R1, --R2, --longreads,
    --longread-signature, --barcode_flavor, --species, --puck,  --puck_barcode_file_id,
    --puck_barcode_file, --investigator, --experiment, --sequencing_date,
    --run_mode arguments.

    :param species_required: if true, the --species argument will be required
        during parsing.
    :param reads_required: if true, --R1, --R2, and --longreads arguments will be
        required during parsing.
    """
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)
    logger.debug("get_data_parser() called")
    parser.add_argument(
        "--R1",
        type=str,
        help=".fastq.gz file path to R1 reads",
        nargs="+",
        required=reads_required,
    )

    parser.add_argument(
        "--R2",
        type=str,
        help=".fastq.gz file path to R2 reads",
        nargs="+",
        required=reads_required,
    )

    parser.add_argument(
        "--reads",
        type=str,
        default="None",
        help="path to a CSV file listing reads for the sample (can be used for pooling reads and also for aggregating bulk samples into one big analysis with different cell barcodes for each sample)",
        required=reads_required,
    )

    parser.add_argument(
        "--dge",
        type=str,
        help="Path to dge matrix. spacemake can also handle already processed"
        + " digital expression data",
        required=reads_required,
    )

    parser.add_argument(
        "--longreads",
        type=str,
        help="fastq(.gz)|fq(.gz)|bam file path to pacbio long reads for library debugging",
        required=reads_required,
    )

    parser.add_argument(
        "--longread-signature",
        type=str,
        help="identify the expected longread signature (see longread.yaml)",
    )

    return parser


def get_set_remove_variable_subparsers(
    parent_parser, variable_name, func, prepend="", allow_multiple=False
):
    """get_set_remove_variable_subparsers.

    :param parent_parser:
    :param variable_name:
    :param func:
    :param prepend:
    :param allow_multiple:
    """
    if allow_multiple:
        nargs = "+"
    else:
        nargs = None

    logger.debug(f"get_set_remove_variable_subparsers(varname={variable_name}) called")

    def get_action_parser(action):
        """Create and return a parser for a specific action.

        :param action: The action (e.g., "set" or "remove")
        """
        # Convert `action` and `variable_name` to kebab-case for subcommand
        kebab_action_variable = f"{action}-{variable_name.replace('_', '-')}"
        snake_action_variable = f"{action}_{variable_name}"

        # Add the parser for the action
        action_parser = parent_parser.add_parser(
            kebab_action_variable,
            parents=[get_project_sample_parser(allow_multiple=True)],
            description=f"{action} {variable_name} for several projects/samples",
            help=f"{action} {variable_name} for several projects/samples",
        )

        # Add the kebab-case and legacy snake_case arguments
        action_parser.add_argument(
            f"--{prepend}{variable_name.replace('_', '-')}",
            f"--{prepend}{variable_name}",  # Legacy argument
            type=str,
            required=True,
            nargs=nargs,
            help=f"Specify the {variable_name} to {action}",
            dest=variable_name,  # Unified destination
        )

        action_parser.set_defaults(func=func, variable=variable_name, action=action)

        return action_parser

    set_parser = get_action_parser("set")

    if allow_multiple:
        set_parser.add_argument("--keep-old", action="store_true")

        remove_parser = get_action_parser("remove")


def get_action_sample_parser(parent_parser, action, func):
    """get_action_sample_parser.

    :param parent_parser:
    :param action:
    :param func:
    """
    logger.debug(f"get_action_sample_parser(action={action}) called")

    if action not in ["add", "update", "delete", "merge"]:
        raise ValueError(f"Invalid action: {action}")

    if action == "merge":
        parser_name = "merge-samples"
        msg = "merge samples"
        parents = [
            get_project_sample_parser(
                prepend="merged-",
                help_extra="of the newly created merged sample"
            ),
            get_project_sample_parser(
                allow_multiple=True,
                help_extra="of the samples to be merged"
            ),
        ]
    else:
        parser_name = f"{action}-sample"
        msg = f"{action} a sample"
        parents = [get_project_sample_parser()]

    # Add additional arguments based on action
    if action == "add":
        # add arguments for species, run_mode, barcode_flavor and puck
        parents.append(
            get_sample_main_variables_parser(
                defaults=True,
                species_required=True,
            )
        )
        # add arguments for R1/R1, dge, longread
        parents.append(get_data_parser())
        # add arguments for extra sample info
        parents.append(get_sample_extra_info_parser())
    elif action == "update":
        # add main variables parser
        parents.append(get_sample_main_variables_parser())
        # add arguments for R1/R1, dge, longread
        parents.append(get_data_parser())
        # add possibility to add extra info
        parents.append(get_sample_extra_info_parser())
    elif action == "merge":
        # add main variables parser
        # when merging, only let the user overwrite puck and run_mode
        parents.append(
            get_sample_main_variables_parser(
                main_variables=["run_mode", "puck"],
            )
        )
        
        # add possibility to add extra info
        parents.append(get_sample_extra_info_parser())

    # Add the parser with kebab-case name and backward-compatible snake_case alias
    sample_parser = parent_parser.add_parser(
        parser_name,
        description=msg,
        help=msg,
        parents=parents,
        aliases=[f"{action}_sample"],
    )

    sample_parser.set_defaults(func=func, action=action)


def setup_project_parser(parent_parser_subparsers):
    """setup_project_parser.

    :param attach_to:
    """
    parser = parent_parser_subparsers.add_parser(
        "projects",
        help="manage projects and samples",
        description="Using one of the subcommands specified below, it is possible to"
        + " add/update/remove projects and their settings",
    )
    subparsers = parser.add_subparsers()

    help_desc = {
        "add_sample_sheet": "add projects and samples from Illumina sample sheet",
        "add_samples_from_yaml": "add several samples at once from a .yaml file",
        "merge_samples": "merge several samples into one. "
        + "samples need to have the same species",
        "list": "list several project(s) and sample(s) and their settings",
    }

    # ADD/UPDATE/DELETE/MERGE sample(s)
    # get parsers for adding, deleting and updating a sample
    for action in ["add", "update", "delete"]:
        get_action_sample_parser(
            subparsers,
            action,
            # attach the function to the parser
            func=add_update_delete_sample_cmdline,
        )

    # get parser for merging
    get_action_sample_parser(subparsers, "merge", func=merge_samples_cmdline)

    # LIST PROJECTS
    # always show these variables
    always_show = [
        "species",
        "barcode_flavor",
        "run_mode",
        "map_strategy",
        "puck",
        "puck_barcode_file_id",
    ]
    import spacemake.project_df as spdf

    remaining_options = [
        x
        for x in spdf.ProjectDF.project_df_default_values.keys()
        if x not in always_show
    ]
    list_projects = subparsers.add_parser(
        "list",
        description=help_desc["list"],
        help=help_desc["list"],
        parents=[
            get_project_sample_parser(
                allow_multiple=True,
                help_extra="subset the data. If none provided, all will be displayed",
            )
        ],
    )
    list_projects.add_argument(
        "--variables",
        help="which extra variables to display per sample? "
        + f"{always_show} will always be shown.",
        choices=remaining_options,
        default=[],
        nargs="*",
    )
    list_projects.set_defaults(func=list_projects_cmdline, always_show=always_show)

    # ADD SAMPLE SHEET
    sample_add_sample_sheet = subparsers.add_parser(
        "add-sample-sheet",
        description=help_desc["add_sample_sheet"],
        help=help_desc["add_sample_sheet"],
        parents=[get_add_sample_sheet_parser()],
    )
    sample_add_sample_sheet_legacy = subparsers.add_parser(
        "add_sample_sheet",
        parents=[get_add_sample_sheet_parser()],
    )
    sample_add_sample_sheet.set_defaults(func=add_sample_sheet_cmdline)
    sample_add_sample_sheet_legacy.set_defaults(func=add_sample_sheet_cmdline)

    # ADD SAMPLES FROM YAML
    sample_add_samples_yaml = subparsers.add_parser(
        "add-samples-from-yaml",
        description=help_desc["add_samples_from_yaml"],
        help=help_desc["add_samples_from_yaml"],
    )
    sample_add_samples_yaml_legacy = subparsers.add_parser("add_samples_from_yaml",)
    sample_add_samples_yaml.add_argument(
        "--samples_yaml",
        type=str,
        required=True,
        help="path to the .yaml file containing sample info",
    )
    sample_add_samples_yaml_legacy.add_argument(
        "--samples_yaml",
        type=str,
        required=True,
        help="path to the .yaml file containing sample info",
    )
    sample_add_samples_yaml.set_defaults(func=add_samples_from_yaml_cmdline)
    sample_add_samples_yaml_legacy.set_defaults(func=add_samples_from_yaml_cmdline)

    # get set/remove parser for each main variable
    # this will add parser for:
    # setting/removing run_mode
    # setting species
    # setting puck
    # setting barcode_flavor
    # NOTE: the ConfigFile.main_variable_sg2type keys determine
    # which commandline args will be defined for the parser.
    # this is a little un-intuitive...
    # TODO: cleaner factory functions for the commandline-parsers
    for main_var_sg, main_var_type in ConfigFile.main_variable_sg2type.items():
        allow_multiple = False
        if isinstance(main_var_type, str) and main_var_type.endswith("_list"):
            allow_multiple = True

        get_set_remove_variable_subparsers(
            subparsers,
            variable_name=main_var_sg,
            func=set_remove_variable_cmdline,
            allow_multiple=allow_multiple,
        )
    return parser


def get_run_parser():
    """get_run_parser."""
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)

    parser.add_argument(
        "--cores", default=1, type=int, help="number of cores to be used in total"
    )
    parser.add_argument(
        "--dryrun",
        "-n",
        action="store_true",
        help="invokes a dry snakemake run, printing only commands",
    )
    parser.add_argument(
        "--debug",
        default="",
        help=f"comma-separated list of logging-domains for which you want DEBUG output",
    )
    parser.add_argument(
        "--rerun-incomplete",
        "--ri",
        action="store_true",
        help="forces snakemake to rerun incompletely generated files",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="if a job fails, keep executing independent jobs",
    )
    parser.add_argument(
        "--printshellcmds",
        "-p",
        action="store_true",
        help="print shell commands for each rule, if exist",
    )
    parser.add_argument(
        "--touch",
        "-t",
        action="store_true",
        help="rather than running the rules, just touch each file",
    )
    
    parser.add_argument(
        "--with-fastqc",
        "-wfqc",
        action="store_true",
        help="Run also fastqc as part of the spacemake run",
        dest="with_fastqc",
    )
    parser.add_argument(
        "--with_fastqc",
        action="store_true",
        help=argparse.SUPPRESS,
        dest="with_fastqc",
    )

    return parser


def setup_init_parser(parent_parser_subparsers):
    """setup_init_parser.

    :param parent_parser_subparsers:
    """
    parser_init = parent_parser_subparsers.add_parser(
        "init",
        help="initialise spacemake: create config files, download genomes and annotations",
    )

    parser_init.add_argument(
        "--root-dir",
        default="",
        help="where to output the results of the spacemake run. defaults to .",
        dest="root_dir",
    )
    parser_init.add_argument(
        "--root_dir",
        default="",
        help=argparse.SUPPRESS,
        dest="root_dir",
    )
    
    
    parser_init.add_argument(
        "--temp-dir",
        default="/tmp",
        help="temporary directory used when running spacemake. defaults to /tmp",
        dest="temp_dir",
    )
    parser_init.add_argument(
        "--temp_dir",
        default="/tmp",
        help=argparse.SUPPRESS,
        dest="temp_dir",
    )
    
    parser_init.add_argument(
        "--download-species",
        default=False,
        help="if set, upon initialisation, spacemake will download the mouse and human genome and index",
        action="store_true",
        dest="download_species",
    )
    parser_init.add_argument(
        "--download_species",
        default=False,
        help=argparse.SUPPRESS,
        action="store_true",
        dest="download_species",
    )
    
    parser_init.add_argument(
        "--dropseq-tools",
        help="absolute path to dropseq_tools directory",
        required=True,
        dest="dropseq_tools",
    )
    parser_init.add_argument(
        "--dropseq_tools",
        help=argparse.SUPPRESS,
        required=False,
        dest="dropseq_tools",
    )

    parser_init.set_defaults(func=spacemake_init)

    return parser_init


def setup_run_parser(parent_parser_subparsers):
    """setup_run_parser.

    :param parent_parser_subparsers:
    """
    parser_run = parent_parser_subparsers.add_parser(
        "run", help="run spacemake", parents=[get_run_parser()]
    )
    parser_run.set_defaults(func=spacemake_run)

    # create subparsers
    parser_run_subparsers = parser_run.add_subparsers()

    downsampling_parser = parser_run_subparsers.add_parser(
        "downsample",
        help="run downsampling analysis for a list of projects/samples",
        parents=[get_project_sample_parser(allow_multiple=True), get_run_parser()],
    )
    downsampling_parser.set_defaults(downsample=True, func=spacemake_run)

    # parser_novosparc = novosparc_spacemake_parser(parser_run_subparsers)
    # parser_novosparc.set_defaults(novosparc_reconstruct=True,
    #    func=lambda args: spacemake_run(pdf, args))

    return parser_run


#####################################################
# actual command-line functions, used as call-backs #
#####################################################


@message_aggregation(logger_name)
def spacemake_init(args):
    """spacemake_init.

    :param args:
    """
    import snakemake
    import yaml
    from spacemake.snakemake.variables import config_path

    if os.path.isfile(config_path):
        msg = "spacemake has already been initiated in this directory.\n"
        msg += "use other commands to run and analyse samples."
        logger.info(msg)
        return 0
    initial_config = os.path.join(os.path.dirname(__file__), "data/config/config.yaml")

    # initialise config file
    cf = ConfigFile.from_yaml(initial_config)
    # update the file path
    cf.set_file_path(config_path)

    # add variables from args
    cf.variables["root_dir"] = args["root_dir"]
    cf.variables["temp_dir"] = args["temp_dir"]

    cf.variables["external_bin"] = {"dropseq_tools": args["dropseq_tools"]}

    cf.variables["microscopy_out"] = args.get("microscopy_out", "")

    if args["download_species"]:
        species_data_config_file = os.path.join(
            os.path.dirname(__file__), "data/config/species_data_url.yaml"
        )
        species_data_config = yaml.load(
            open(species_data_config_file), Loader=yaml.FullLoader
        )

        # add keys as species to config
        species = list(species_data_config.keys())
        snakemake_config = {"root_dir": ""}
        snakemake_config["species"] = species

        # define the pattern
        snakemake_config["annotation_file_pattern"] = (
            "species_data/{species}/{species}_{data_type}.gtf"
        )
        snakemake_config["genome_file_pattern"] = (
            "species_data/{species}/{species}_{data_type}.fa"
        )

        # the to be saved file paths
        species_info = {}
        species_files_exist = []

        for sp in species:
            if not sp in species_info.keys():
                species_info[sp] = {}
            for data_type in ["genome", "annotation"]:
                # save download url to be passed to snakemake
                snakemake_config[sp + "_" + data_type + "_url"] = species_data_config[
                    sp
                ][data_type]
                # save file path of the result
                species_file = snakemake_config[data_type + "_file_pattern"].format(
                    species=sp, data_type=data_type
                )
                species_info[sp][data_type] = species_file

                # add bool if file exists
                species_files_exist.append(os.path.isfile(species_file))

        # if any of the files are missing
        if not all(species_files_exist):
            # get the snakefile for downloading the .gtf and .fa files
            snakefile = os.path.join(
                os.path.dirname(__file__), "snakemake/species_init.smk"
            )
            # run snakemake: download species data and place them in the right place
            snakemake.snakemake(snakefile, cores=1, config=snakemake_config)

        for key, value in species_info.items():
            cf.add_variable(
                "species",
                key,
                reference="genome",
                sequence=value["genome"],
                annotation=value["annotation"],
            )

    # copy visium_puck_barcode_file
    dest_visium_path = "puck_data/visium_barcode_positions.csv"
    logger.info(f"Moving visium barcodes to {dest_visium_path}")
    os.makedirs(os.path.dirname(dest_visium_path), exist_ok=True)
    from shutil import copyfile

    copyfile(
        os.path.join(os.path.dirname(__file__), "data/visium_barcode_positions.csv"),
        dest_visium_path,
    )

    # copy openst_coordinate_system
    dest_puck_collection_path = "puck_data/openst_coordinate_system.csv"
    logger.info(
        f"Moving puck collection coordinate system to {dest_puck_collection_path}"
    )
    os.makedirs(os.path.dirname(dest_puck_collection_path), exist_ok=True)
    copyfile(
        os.path.join(
            os.path.dirname(__file__),
            "data/puck_collection/openst_coordinate_system.csv",
        ),
        dest_puck_collection_path,
    )

    # save
    cf.dump()


@message_aggregation(logger_name)
def spacemake_run(args):
    """spacemake_run.

    :param args:
    """
    from spacemake.errors import SpacemakeError
    import snakemake

    import spacemake.snakemake.variables as var

    if not os.path.isfile(var.config_path):
        msg = "spacemake has not been initalised yet.\n"
        msg += "please run `spacemake init` to start a new project"

        raise SpacemakeError(msg)

    from spacemake.project_df import get_global_ProjectDF

    root = logging.getLogger("spacemake")
    # TODO: Why is args a dictionary and not argparse.Namespace ?
    if args["debug"]:
        # activate cmdline requested debug output for specific domains (comma-separated)
        for logger_name in args["debug"].split(","):
            if logger_name:
                root.info(f"setting domain {logger_name} to DEBUG")
                logging.getLogger(logger_name.replace("root", "")).setLevel(
                    logging.DEBUG
                )

    pdf = get_global_ProjectDF()
    samples = []
    projects = []
    targets = ["run_analysis"]
    with_fastqc = args.get("with_fastqc", False)

    downsample = args.get("downsample", False)
    novosparc_reconstruct = args.get("novosparc_reconstruct", False)

    if downsample:
        targets = ["downsample"]
        samples = args.get("sample_id_list", [])
        projects = args.get("project_id_list", [])

    if novosparc_reconstruct:
        targets = ["novosparc"]
        import spacemake.smk as smk

        novosparc_variables = smk.get_novosparc_variables(pdf, args)
    else:
        novosparc_variables = {}

    config_variables = {
        "project_df": pdf.file_path,
        "samples": samples,
        "projects": projects,
        "with_fastqc": with_fastqc,
        "pwd": os.getcwd(),
        "log_debug": args["debug"],
    }

    # join config_variables and novosparc_variables
    # to flatten the dictionary
    config_variables = {**config_variables, **novosparc_variables}

    # get the snakefile
    snakefile = os.path.join(os.path.dirname(__file__), "snakemake/main.smk")
    # run snakemake
    preprocess_finished = snakemake.snakemake(
        snakefile,
        configfiles=[var.config_path],
        cores=args["cores"],
        dryrun=args["dryrun"],
        targets=["get_stats_prealigned_barcodes", "unload_genome_flag"],
        touch=args["touch"],
        force_incomplete=args["rerun_incomplete"],
        keepgoing=args["keep_going"],
        printshellcmds=args["printshellcmds"],
        config=config_variables,
    )

    if preprocess_finished is False:
        raise SpacemakeError("an error occurred while snakemake() ran")

    # update valid pucks (above threshold) before continuing to downstream
    # this performs counting of matching barcodes after alignment
    pdf.update_project_df_barcode_matches(prealigned=True)
    pdf.consolidate_pucks_merged_samples()
    pdf.dump()

    # whitelisting of barcodes
    preprocess_finished = snakemake.snakemake(
        snakefile,
        configfiles=[var.config_path],
        cores=args["cores"],
        dryrun=args["dryrun"],
        targets=["get_whitelist_barcodes"],
        touch=args["touch"],
        force_incomplete=args["rerun_incomplete"],
        keepgoing=args["keep_going"],
        printshellcmds=args["printshellcmds"],
        config=config_variables,
    )

    pdf.update_project_df_barcode_matches()
    pdf.dump()

    # run snakemake quantification and reports
    analysis_finished = snakemake.snakemake(
        snakefile,
        configfiles=[var.config_path],
        cores=args["cores"],
        dryrun=args["dryrun"],
        targets=targets,
        touch=args["touch"],
        force_incomplete=args["rerun_incomplete"],
        keepgoing=args["keep_going"],
        printshellcmds=args["printshellcmds"],
        config=config_variables,
    )

    if analysis_finished is False:
        raise SpacemakeError("an error occurred while snakemake() ran")

    # at the very end dump the project_data_frame
    pdf.dump()


@message_aggregation(logger_name)
def add_update_delete_sample_cmdline(args):
    """add_update_delete_sample_cmdline.

    :param args:
    """
    action = args["action"]
    from spacemake.project_df import get_global_ProjectDF

    pdf = get_global_ProjectDF()
    # remove the action from args
    del args["action"]

    if action == "add" or action == "update":
        func = lambda **kwargs: pdf.add_update_sample(action=action, **kwargs)
    elif action == "delete":
        func = pdf.delete_sample

    sample = func(**args)
    pdf.dump()


@message_aggregation(logger_name)
def add_sample_sheet_cmdline(args):
    """add_sample_sheet_cmdline.

    :param args:
    """
    from spacemake.project_df import get_global_ProjectDF

    pdf = get_global_ProjectDF()

    pdf.add_sample_sheet(args["sample_sheet"], args["basecalls_dir"])

    pdf.dump()


@message_aggregation(logger_name)
def add_samples_from_yaml_cmdline(args):
    """add_samples_from_yaml_cmdline.

    :param args:
    """
    from spacemake.project_df import get_global_ProjectDF

    pdf = get_global_ProjectDF()
    pdf.add_samples_from_yaml(args["samples_yaml"])

    pdf.dump()


@message_aggregation(logger_name)
def set_remove_variable_cmdline(args):
    """set_remove_variable_cmdline.

    :param args:
    """
    variable_name = args["variable"]
    action = args["action"]

    variable_key = args[variable_name]
    from spacemake.project_df import get_global_ProjectDF

    pdf = get_global_ProjectDF()

    pdf.set_remove_variable(
        variable_name=variable_name,
        variable_key=variable_key,
        action=action,
        project_id_list=args["project_id_list"],
        sample_id_list=args["sample_id_list"],
        keep_old=args.get("keep_old", False),
    )

    pdf.dump()


@message_aggregation(logger_name)
def merge_samples_cmdline(args):
    """merge_samples_cmdline.

    :param args:
    """
    from spacemake.project_df import get_global_ProjectDF

    pdf = get_global_ProjectDF()
    pdf.merge_samples(**args)

    pdf.dump()


@message_aggregation(logger_name)
def list_projects_cmdline(args):
    """list_projects_cmdline.

    :param args:
    """
    from spacemake.project_df import get_global_ProjectDF

    pdf = get_global_ProjectDF()

    projects = args["project_id_list"]
    samples = args["sample_id_list"]
    variables = args["always_show"] + args["variables"]

    df = pdf.df
    logger = logging.getLogger(logger_name)

    if projects != [] or samples != []:
        df = df.query("project_id in @projects or sample_id in @samples")
        logger.inof(f"listing projects: {projects} and samples: {samples}")
    else:
        logger.info("listing all projects and samples")

    logger.info(f"variables used: {variables}")

    # print the table
    logger.info(df.loc[:, variables].__str__())


def make_main_parser():
    #################
    # DEFINE PARSER #
    #################
    # we add main parser to a dictionary
    # so that we can later call the help function
    # based on the sub-command. this is to print the
    # -h (help) functionality if no parameters are provided,
    # rather than printing an error.
    # import spacemake.smk as smk

    parser_main = argparse.ArgumentParser(
        allow_abbrev=False,
        description="spacemake: bioinformatic pipeline for processing and analysis of spatial-transcriptomics data",
    )

    parser_main.add_argument("--version", action="store_true")

    parser_main_subparsers = parser_main.add_subparsers(
        help="sub-command help", dest="subcommand"
    )

    parser_run = None
    parser_projects = None
    parser_config = None
    parser_init = None
    parser_spatial = None

    ##################
    # SPACEMAKE INIT #
    ##################
    parser_init = setup_init_parser(parser_main_subparsers)

    ####################
    # SPACEMAKE CONFIG #
    ####################
    ## spacemake_config args
    from spacemake.config import setup_config_parser
    from spacemake.snakemake.variables import config_path

    if os.path.isfile(config_path):
        # save config file
        parser_config = setup_config_parser(parser_main_subparsers)

        ############################
        # SPACEMAKE PROJECT/SAMPLE #
        ############################
        from spacemake.cmdline import setup_project_parser

        parser_projects = setup_project_parser(parser_main_subparsers)

        #################
        # SPACEMAKE RUN #
        #################
        parser_run = setup_run_parser(parser_main_subparsers)

        #####################
        # SPACEMAKE SPATIAL #
        #####################
        from spacemake.spatial.cmdline import setup_spatial_parser

        parser_spatial = setup_spatial_parser(parser_main_subparsers)

    parser_dict = {
        "init": parser_init,
        "config": parser_config,
        "projects": parser_projects,
        "run": parser_run,
        "main": parser_main,
        "spatial": parser_spatial,
    }

    return parser_main, parser_dict


def cmdline():
    # import importlib.metadata
    """cmdline."""
    parser_main, parser_dict = make_main_parser()
    args = parser_main.parse_args()

    if args.version and args.subcommand is None:
        import spacemake.contrib

        print(spacemake.contrib.__version__)
        return 0
    else:
        del args.version

    # get the function to be run
    if "func" in args:
        func = args.func
    # else print help
    else:
        if args.subcommand is not None:
            parser_dict[args.subcommand].print_help()
        else:
            parser_dict["main"].print_help()
        return 0

    # get the args and delete the func key, get only set values
    args = {key: value for key, value in vars(args).items() if value is not None}
    args.pop("func", None)
    # pop also main,
    args.pop("subcommand", None)

    return func(args)


if __name__ == "__main__":
    cmdline()
