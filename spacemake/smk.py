import re
import sys
import os
import snakemake
import argparse
import yaml
import logging
import scanpy as sc
import pandas as pd


from shutil import copyfile
from spacemake.project_df import ProjectDF, get_project_sample_parser
from spacemake.config import ConfigFile
from spacemake.util import message_aggregation
from spacemake.errors import SpacemakeError
from spacemake.preprocess import attach_puck_variables

config_path = "config.yaml"
project_df = "project_df.csv"

logger_name = "spacemake.main"
logger = logging.getLogger(logger_name)

class Spacemake:
    def __init__(self, root):
        """__init__.

        :param root:
        """
        self.root = root
        self.config = ConfigFile.from_yaml(f'{root}/config.yaml')
        self.project_df = ProjectDF(f'{root}/project_df.csv', config=self.config)
    
    def load_processed_adata(self,
        project_id,
        sample_id,
        run_mode_name,
        umi_cutoff 
    ):
        """load_processed_adata.

        :param project_id:
        :param sample_id:
        :param run_mode_name:
        :param umi_cutoff:
        """
        self.project_df.assert_run_mode(project_id, sample_id, run_mode_name)
        run_mode = self.config.get_run_mode(run_mode_name)

        if not int(umi_cutoff) in [int(uc)
            for uc in run_mode.variables['umi_cutoff']]:
            raise SpacemakeError(f'run_mode={run_mode} has no ' + 
                f'umi_cutoff={umi_cutoff}')

        adata_raw = self.load_raw_spatial_adata(
            project_id = project_id,
            sample_id = sample_id,
            run_mode_name = run_mode_name)

        adata = sc.read(f'{self.root}/projects/{project_id}/processed_data/{sample_id}/'+
                f'illumina/complete_data/automated_analysis/{run_mode_name}/' +
                f'umi_cutoff_{umi_cutoff}/results.h5ad')
        
        if 'run_mode_variables' not in adata.uns.keys():
            adata.uns['run_mode_variables'] = run_mode.variables
        if 'puck_variables' not in adata.uns.keys():
            adata.uns['puck_variables'] = adata_raw.uns['puck_variables']
        
        return adata
    
    def load_raw_spatial_adata(
        self,
        project_id,
        sample_id,
        run_mode_name
    ):
        """load_raw_spatial_adata.

        :param project_id:
        :param sample_id:
        :param run_mode_name:
        """
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

        adata = sc.read(f'{self.root}/projects/{project_id}/processed_data/{sample_id}/'+
                f'illumina/complete_data/dge/dge{dge_type}{dge_cleaned}'+
                f'{polyA_adapter_trimmed}{mm_included}.spatial_beads.h5ad')
        
        if 'puck_variables' not in adata.uns.keys():
            adata = attach_puck_variables(adata,
                puck_variables = self.project_df.get_puck_variables(
                    project_id = project_id,
                    sample_id = sample_id
                )
            )

        if 'run_mode_variables' not in adata.uns.keys():
            adata.uns['run_mode_variables'] = run_mode.variables

        return adata

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

def setup_run_parser(pdf, parent_parser):
    parser_run = parent_parser.add_parser(
        "run", help="run spacemake", parents=[get_run_parser()]
    )
    parser_run.set_defaults(func=lambda args: spacemake_run(pdf, args))

    # create subparsers
    parser_run_subparsers = parser_run.add_subparsers()

    downsampling_parser = parser_run_subparsers.add_parser(
        "downsample",
        help="run downsampling analysis for a list of projects/samples",
        parents=[get_project_sample_parser(allow_multiple=True), get_run_parser()],
    )
    downsampling_parser.set_defaults(downsample=True,
        func=lambda args: spacemake_run(pdf, args))

    parser_novosparc = parser_run_subparsers.add_parser(
        'novosparc',
        help='reconstruct the 2-d tissue with novosparc',
        parents=[
            get_run_parser(),
            get_project_sample_parser(
                    allow_multiple = False,
                    help_extra = ' of sample to reconstruct',
                 ),
        ],
    )
        
    parser_novosparc.add_argument(
        '--run_mode', type=str, help='the run_mode to be used for this sample ' +
        'for reconstruction. If left empty, the first run_mode for this ' +
        'sample will be used', required=False,
    )

    parser_novosparc.add_argument(
        '--umi_cutoff', type=int, help='umi_cutoff to be used for this sample ' +
        'for reconstruction. If left empty, the smallest umi_cutoff of a given ' +
        'run_mode will be used', required=False,
    )

    parser_novosparc.add_argument(
        '--reference_project_id', type=str, help='project_id of the reference atlas.' +
        'Has to be spatial, otherwise error will be raised', required=False,
    )

    parser_novosparc.add_argument(
        '--reference_sample_id', type=str, help='sample_id of the reference atlas.' +
        'Has to be spatial, otherwise error will be raised', required=False,
    )

    parser_novosparc.add_argument(
        '--reference_run_mode', type=str, help='the run_mode to be used for the ' +
        'reference sample. If empty, the first run_mode for the reference will be used',
        required=False,
    )

    parser_novosparc.add_argument(
        '--reference_umi_cutoff', type=str, help='the umi_cutoff to be used for this ' +
        'reference sample. If empty, the smallest umi_cutoff of the given run_mode ' +
        'will be used.', required=False,
    )

    parser_novosparc.set_defaults(novosparc_reconstruct=True,
        func=lambda args: spacemake_run(pdf, args))
    
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

def get_novosparc_variables(pdf, args):
    # assert that sample exists
    pdf.assert_sample(args['project_id'], args['sample_id'])

    def populate_variables_from_args(pdf, args, arg_prefix=''):
        # get sample info
        sample_info = pdf.get_sample_info(
            project_id=args[f'{arg_prefix}project_id'],
            sample_id=args[f'{arg_prefix}sample_id']
        )

        # populate return dictionary
        ret = {
            f'{arg_prefix}project_id': args[f'{arg_prefix}project_id'],
            f'{arg_prefix}sample_id': args[f'{arg_prefix}sample_id'],
        }

        # get run mode
        if f'{arg_prefix}run_mode' in args:
            ret[f'{arg_prefix}run_mode'] = args[f'{arg_prefix}run_mode']
        else:
            run_mode_name = sample_info['run_mode'][0]
            ret[f'{arg_prefix}run_mode'] = run_mode_name
            logger.info(f"No run_mode provided, using {run_mode_name}")

        run_mode = pdf.config.get_run_mode(ret[f'{arg_prefix}run_mode'])

        if f'{arg_prefix}umi_cutoff' not in args:
            umi_cutoff = run_mode.variables['umi_cutoff'][0]
            ret[f'{arg_prefix}umi_cutoff'] = umi_cutoff
            logger.info(f"No umi_cutoff provided, using {umi_cutoff}")
        else:
            ret[f'{arg_prefix}umi_cutoff'] = args[f'{arg_prefix}umi_cutoff']

        return ret
    
    ret = populate_variables_from_args(pdf, args)

    if 'reference_project_id' not in args or 'reference_sample_id' not in args:
        logger.info('No reference_project_id or reference_sample_id provided,' +
            ' running novosparc de-novo...')
        ret['reference_project_id'] = ''
        ret['reference_sample_id'] = ''
        ret['reference_umi_cutoff'] = ''
        ret['reference_run_mode'] = ''
    else:
        pdf.assert_sample(args['reference_project_id'],
            args['reference_sample_id'])

        logger.info("Using (project_id, sample_id)="+
            f"({args['reference_project_id']}, {args['reference_sample_id']})" +
            " reference, running novosparc with reference...")

        novosparc_ret = populate_variables_from_args(pdf, args,
            arg_prefix='reference_')

        ret = {**ret, **novosparc_ret}

    return ret

@message_aggregation(logger_name)
def spacemake_run(pdf, args):
    if not os.path.isfile(config_path):
        msg = "spacemake has not been initalised yet.\n"
        msg += "please run `spacemake init` to start a new project"

        raise SpacemakeError(msg)

    samples = []
    projects = []
    targets = ['all']
    with_fastqc = args.get("with_fastqc", False)

    downsample = args.get("downsample", False)
    novosparc_reconstruct = args.get('novosparc_reconstruct', False)

    if downsample:
        targets = ["downsample"]
        samples = args.get("sample_id_list", [])
        projects = args.get("project_id_list", [])

    if novosparc_reconstruct:
        targets = ["novosparc"]
        novosparc_variables = get_novosparc_variables(pdf, args)
    else:
        novosparc_variables = {}

    config_variables = {
        "project_df": pdf.file_path,
        "samples": samples,
        "projects": projects,
        "with_fastqc": with_fastqc,
    }

    # join config_variables and novosparc_variables
    # to flatten the directory
    config_variables = {**config_variables, **novosparc_variables}

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
        config=config_variables,
    )

    if run_successful is False:
        raise SpacemakeError("an error occurred while snakemake() ran")

    # at the very end dump the project_data_frame
    pdf.dump()

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
    parser_run = setup_run_parser(pdf, subparsers)

#####################
# SPACEMAKE SPATIAL #
#####################
from .spatial.cmdline import setup_spatial_parser

if os.path.isfile(config_path):
    spmk = Spacemake('.')
    parser_spatial = setup_spatial_parser(spmk, subparsers)

def cmdline():
    args = main_parser.parse_args()

    parser_dict = {
        "init": parser_init,
        "config": parser_config,
        "projects": parser_projects,
        "run": parser_run,
        "main": main_parser,
        "spatial": parser_spatial
    }

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

    func(args)

if __name__ == "__main__":
    cmdline()
