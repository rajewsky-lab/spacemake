import re
import sys
import os
import snakemake
import argparse
import yaml
import logging

from shutil import copyfile
from spacemake.project_df import ProjectDF, get_project_sample_parser
from spacemake.config import ConfigFile
from spacemake.util import message_aggregation
from spacemake.errors import SpacemakeError

config_path = "config.yaml"
project_df = "project_df.csv"

logger_name = "spacemake.main"
logger = logging.getLogger(logger_name)


def get_run_parser():
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
        "--with_fastqc",
        "-wfqc",
        action="store_true",
        help="Run also fastqc as part of the spacemake run",
    )

    return parser


def setup_init_parser(parent_parser):
    parser_init = parent_parser.add_parser(
        "init",
        help="initialise spacemake: create config files, download genomes and annotations",
    )
    parser_init.add_argument(
        "--root_dir",
        default="",
        help="where to output the results of the spacemake run. defaults to .",
    )
    parser_init.add_argument(
        "--temp_dir",
        default="/tmp",
        help="temporary directory used when running spacemake. defaults to /tmp",
    )
    parser_init.add_argument(
        "--download_species",
        default=False,
        help="if set, upon initialisation, spacemake will download the mouse and human genome and index",
        action="store_true",
    )
    parser_init.add_argument(
        "--dropseq_tools",
        help="absolute path to dropseq_tools directory",
        required=True,
    )
    parser_init.set_defaults(func=spacemake_init)

    return parser_init


def setup_run_parser(parent_parser):
    parser_run = parent_parser.add_parser(
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

    return parser_run


@message_aggregation(logger_name)
def spacemake_init(args):
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

    cf.variables['external_bin'] = {
        'dropseq_tools' : args['dropseq_tools']
    }

    cf.variables['microscopy_out'] = args.get('microscopy_out', '')

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
        snakemake_config[
            "annotation_file_pattern"
        ] = "species_data/{species}/{species}_{data_type}.gtf"
        snakemake_config[
            "genome_file_pattern"
        ] = "species_data/{species}/{species}_{data_type}.fa"

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
                "species", key, genome=value["genome"], annotation=value["annotation"]
            )

    # copy visium_puck_barcode_file
    dest_visium_path = "puck_data/visium_barcode_positions.csv"
    logger.info(f"Moving visium barcodes to {dest_visium_path}")
    os.makedirs(os.path.dirname(dest_visium_path), exist_ok=True)
    copyfile(
        os.path.join(os.path.dirname(__file__), "data/visium_barcode_positions.csv"),
        dest_visium_path,
    )

    # save
    cf.dump()


@message_aggregation(logger_name)
def spacemake_run(args):
    if not os.path.isfile(config_path):
        msg = "spacemake has not been initalised yet.\n"
        msg += "please run `spacemake init` to start a new project"

        raise SpacemakeError(msg)

    samples = []
    projects = []
    targets = None
    with_fastqc = args.get("with_fastqc", False)

    downsample = args.get("downsample", False)

    if downsample:
        targets = ["downsample"]
        samples = args.get("sample_id_list", [])
        projects = args.get("project_id_list", [])

    # get the snakefile
    snakefile = os.path.join(os.path.dirname(__file__), "snakemake/main.smk")
    # run snakemake
    run_successful = snakemake.snakemake(
        snakefile,
        configfiles=[config_path],
        cores=args["cores"],
        dryrun=args["dryrun"],
        targets=targets,
        touch=args["touch"],
        force_incomplete=args["rerun_incomplete"],
        keepgoing=args["keep_going"],
        printshellcmds=args["printshellcmds"],
        config={
            "project_df": project_df,
            "samples": samples,
            "projects": projects,
            "with_fastqc": with_fastqc,
        },
    )

    if run_successful is False:
        raise SpacemakeError("an error occurred while snakemake() ran")


#################
# DEFINE PARSER #
#################
# we add main parser to a dictionary
# so that we can later call the help function
# based on the sub-command. this is to print the
# -h (help) functionality if no parameters are provided,
# rather than printing an error.

main_parser = argparse.ArgumentParser(
    allow_abbrev=False,
    description="spacemake: bioinformatic pipeline for processing and analysis of spatial-transcriptomics data",
)

subparsers = main_parser.add_subparsers(help="sub-command help", dest="subcommand")

parser_run = None
parser_projects = None
parser_config = None
parser_init = None

##################
# SPACEMAKE INIT #
##################
parser_init = setup_init_parser(subparsers)

####################
# SPACEMAKE CONFIG #
####################
## spacemake_config args
from spacemake.config import setup_config_parser

if os.path.isfile(config_path):
    cf = ConfigFile.from_yaml(config_path)
    # save config file
    parser_config = setup_config_parser(cf, subparsers)

############################
# SPACEMAKE PROJECT/SAMPLE #
############################
from spacemake.project_df import setup_project_parser

if os.path.isfile(config_path):
    pdf = ProjectDF(project_df, cf)
    parser_projects = setup_project_parser(pdf, subparsers)

    #################
    # SPACEMAKE RUN #
    #################
    parser_run = setup_run_parser(subparsers)


def cmdline():
    args = main_parser.parse_args()

    parser_dict = {
        "init": parser_init,
        "config": parser_config,
        "projects": parser_projects,
        "run": parser_run,
        "main": main_parser,
    }

    # get the function to be run
    if "func" in args:
        func = args.func
    # else print help
    else:
        if args.subcommand is not None:
            print(args.subcommand)
            parser_dict[args.subcommand].print_help()
        else:
            parser_dict["main"].print_help()
        return 0

    # get the args and delete the func key, get only set values
    args = {key: value for key, value in vars(args).items() if value is not None}
    args.pop("func", None)
    # pop also main,
    args.pop("subcommand", None)

    func(args)


if __name__ == "__main__":
    cmdline()
