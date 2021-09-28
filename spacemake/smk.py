import re
import sys
import os
import snakemake
import argparse
import yaml

from shutil import copyfile
from spacemake.project_df import ProjectDF
from spacemake.config import ConfigFile

config_path = 'config.yaml'
project_df = 'project_df.csv'

def spacemake_init(args):
    if os.path.isfile(config_path):
        msg = "spacemake has already been initiated in this directory.\n"
        msg += "use other commands to run and analyse samples."
        print(f"{msg}")
        return 0
    initial_config = os.path.join(os.path.dirname(__file__),
            'data/config/config.yaml')

    # initialise config file
    cf = ConfigFile(initial_config)
    # update the file path 
    cf.set_file_path(config_path)
    cf.variables['root_dir'] = args['root_dir']
    cf.variables['temp_dir'] = args['temp_dir']

    if args['download_species']:
        species_data_config_file = os.path.join(os.path.dirname(__file__),
                'data/config/species_data_url.yaml')
        species_data_config = yaml.load(open(species_data_config_file),
                        Loader=yaml.FullLoader)
        
        # add keys as species to config
        species = list(species_data_config.keys())
        snakemake_config = {'root_dir':''}
        snakemake_config['species'] = species

        # define the pattern
        snakemake_config['annotation_file_pattern'] = 'species_data/{species}/{species}_{data_type}.gtf'
        snakemake_config['genome_file_pattern'] = 'species_data/{species}/{species}_{data_type}.fa'

        # the to be saved file paths
        species_info = {}
        species_files_exist = []
        
        for sp in species:
            if not sp in species_info.keys():
                species_info[sp] = {}
            for data_type in ['genome', 'annotation']:
                # save download url to be passed to snakemake
                snakemake_config[sp+'_'+data_type+'_url'] = species_data_config[sp][data_type]
                # save file path of the result
                species_file = snakemake_config[data_type+'_file_pattern']\
                        .format(species = sp, data_type=data_type)
                species_info[sp][data_type] = species_file

                # add bool if file exists
                species_files_exist.append(os.path.isfile(species_file))

        # if any of the files are missing
        if not all(species_files_exist):
            # get the snakefile for downloading the .gtf and .fa files
            snakefile = os.path.join(os.path.dirname(__file__), 'snakemake/species_init.smk')
            # run snakemake: download species data and place them in the right place
            snakemake.snakemake(snakefile, cores = 1, config=snakemake_config)

        for key, value in species_info.items():
            cf.add_variable('species', key, genome=value['genome'],
                    annotation=value['annotation'])

    # copy visium_puck_barcode_file
    dest_visium_path = 'puck_data/visium_barcode_positions.csv'
    os.makedirs(os.path.dirname(dest_visium_path), exist_ok=True)
    copyfile(os.path.join(os.path.dirname(__file__), 'data/visium_barcode_positions.csv'),
        dest_visium_path)

    # save
    cf.dump()


def spacemake_run(args):
    if not os.path.isfile(config_path):
        msg = "spacemake has not been initalised yet.\n"
        msg += "please run `spacemake init` to start a new project"
        return 2

    # get the snakefile
    snakefile = os.path.join(os.path.dirname(__file__), 'snakemake/main.smk')
    # run snakemake
    snakemake.snakemake(snakefile, configfiles=[config_path],
        cores = args['cores'], dryrun=args['dryrun'],
        force_incomplete=args['rerun_incomplete'],
        keepgoing=args['keep_going'],
        config={'project_df': project_df})

#################
# DEFINE PARSER #
#################
# we add main parser to a dictionary
# so that we can later call the help function
# based on the sub-command. this is to print the 
# -h (help) functionality if no parameters are provided,
# rather than printing an error.

main_parser= argparse.ArgumentParser(
        allow_abbrev=False,
        description='spacemake: bioinformatic pipeline for processing and analysis of spatial-transcriptomics data')

subparsers = main_parser.add_subparsers(help='sub-command help', dest='subcommand')

##################
# SPACEMAKE INIT #
##################
parser_init = subparsers.add_parser('init', help = 'initialise spacemake: create config files, download genomes and annotations')
parser_init.add_argument('--root_dir', default='',
    help = 'where to output the results of the spacemake run. defaults to .')
parser_init.add_argument('--temp_dir', default='/tmp',
    help='temporary directory used when running spacemake. defaults to /tmp')
parser_init.add_argument('--download_species', default=False,
    help='if set, upon initialisation, spacemake will download the mouse and human genome and index', action='store_true')
parser_init.add_argument('--dropseq_tools',
    help='absolute path to dropseq_tools directory', required=True)
parser_init.add_argument('--picard_tools',
    help='absolute path to `picard.jar` file', required=True)
parser_init.set_defaults(func=spacemake_init)

#################
# SPACEMAKE RUN #
#################
parser_run = subparsers.add_parser('run', help = 'run spacemake')
parser_run.add_argument('--cores',
    default=1,
    type=int,
    help = 'number of cores to be used in total')
parser_run.add_argument('--dryrun', '-n', action='store_true', help = 'invokes a dry snakemake run, printing only commands')
parser_run.add_argument('--rerun-incomplete', '--ri', action='store_true', help = 'forces snakemake to rerun incompletely generated files')
parser_run.add_argument('--keep-going', action='store_true', help='if a job fails, keep executing independent jobs')
parser_run.set_defaults(func=spacemake_run)
####################
# SPACEMAKE CONFIG #
####################
## spacemake_config args
if os.path.isfile(config_path):
    cf = ConfigFile(config_path)
    # save config file
    cf.dump()
    parser_config = cf.get_subparsers(subparsers)

####################
# SPACEMAKE SAMPLE #
####################
from spacemake.project_df import setup_project_parser

if os.path.isfile(config_path):
    pdf = ProjectDF(project_df, cf)
    parser_projects = setup_project_parser(pdf, subparsers)

def cmdline():
    args = main_parser.parse_args()

    parser_dict = {
        'init': parser_init,
        'config': parser_config,
        'projects': parser_projects,
        'run': parser_run}

    # get the function to be run
    if 'func' in args:
        func = args.func
    # else print help
    else:
        if args.subcommand is not None:
            print(args.subcommand)
            parser_dict[args.subcommand].print_help()
        else:
            parser_dict['main'].print_help()
        return 0
        
    # get the args and delete the func key, get only set values
    args = {key: value for key, value in vars(args).items() if value is not None}
    args.pop('func', None)
    # pop also main, 
    args.pop('subcommand', None)

    func(args)

if __name__ == "__main__":
    cmdline()
