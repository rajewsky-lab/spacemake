import re
import sys
import os
import snakemake
import argparse
import yaml

from shutil import copyfile
from spacemake.cmdline_utils import ProjectDF, ConfigFile

config_path = 'config.yaml'
project_df = 'project_df.csv'

def spacemake_init(args):
#    if os.path.isfile(config_path):
#        msg = "spacemake has already been initiated in this directory.\n"
#        msg += "use other commands to run and analyse samples."
#        print(f"{msg}")
#        return 0
    initial_config = os.path.join(os.path.dirname(__file__),
            'config/config.yaml')
    species_data_config_file = os.path.join(os.path.dirname(__file__),
            'config/species_data_url.yaml')
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
        snakemake.snakemake(snakefile, cores = 1, config=snakemake_config, dryrun=True)

    cf = ConfigFile(initial_config)
    # update the file path 
    cf.set_file_path(config_path)
    # add species info which we just generated
    cf.add_species_info(species_info)
    # save
    cf.dump()


def spacemake_run(args):
    if not os.path.isfile(config_path):
        msg = "spacemake has not been initalised yet.\n"
        msg += "please run `spacemake init` to start a new project"
        return 0

    # get the snakefile
    snakefile = os.path.join(os.path.dirname(__file__), 'snakemake/main.smk')
    # run snakemake
    snakemake.snakemake(snakefile, configfiles=[config_path],
        cores = args['cores'], #dryrun=args.dryrun,
        config={'root_dir': '', 'temp_dir': '/tmp'}, dryrun=True)

#################
# DEFINE PARSER #
#################

parsers = {
    'main': argparse.ArgumentParser(
        description='spacemake: bioinformatic pipeline for processing and analysis of spatial-transcriptomics data')
}

subparsers = parsers['main'].add_subparsers(help='sub-command help', dest='main')

##################
# SPACEMAKE INIT #
##################
parsers['init'] = subparsers.add_parser('init', help = 'initialise spacemake: create config files, download genomes and annotations')
parsers['init'].set_defaults(func=spacemake_init)

## spacemake_run args
parsers['run'] = subparsers.add_parser('run', help = 'run spacemake')
parsers['run'].add_argument('--cores',
    default=1,
    type=int,
    help = 'number of cores to be used in total')
parsers['run'].add_argument('--dryrun', '-n', action='store_true', help = 'invokes a dry snakemake run, printing only commands')
parsers['run'].set_defaults(func=spacemake_run)

####################
# SPACEMAKE CONFIG #
####################
## spacemake_config args
parsers['config'] = ConfigFile.add_config_subparsers(subparsers, config_path)

####################
# SPACEMAKE SAMPLE #
####################
parsers['projects'] = subparsers.add_parser('projects', help ='manage projects and samples')
parser_sample_subparsers = parsers['projects'].add_subparsers(help = 'sample sub-command help')

# ADD SAMPLE SHEET
parser_sample_add_sample_sheet = parser_sample_subparsers.add_parser('add_sample_sheet',
    help = 'add projects and samples from Illumina sample sheet',
    parents=[ProjectDF.get_add_sample_sheet_parser()])
parser_sample_add_sample_sheet.set_defaults(func=ProjectDF.add_sample_sheet_cmdline,
    project_df_file=project_df)

# ADD SAMPLES FROM YAML
parser_sample_add_samples_yaml = parser_sample_subparsers.add_parser('add_samples_from_yaml',
    help = 'add several samples at once from a .yaml file')
parser_sample_add_samples_yaml.add_argument('samples_yaml',
    type=str,
    help='path to the .yaml file containing sample info')
parser_sample_add_samples_yaml.set_defaults(
    func = ProjectDF.add_samples_from_yaml,
    project_df_file=project_df)

# ADD SAMPLE
parser_sample_add = parser_sample_subparsers.add_parser('add_sample',
    help = 'add new sample',
    parents=[ProjectDF.get_add_update_sample_parser(),
             ProjectDF.get_read_species_parser(reads_required=True)])
parser_sample_add.set_defaults(func=ProjectDF.add_sample_cmdline,
    project_df_file=project_df)

# UPDATE SAMPLE
parser_sample_add = parser_sample_subparsers.add_parser('update_sample',
    help = 'update_existing_sample',
    parents=[ProjectDF.get_add_update_sample_parser(),
             ProjectDF.get_read_species_parser(reads_required=False)])
parser_sample_add.set_defaults(func=ProjectDF.update_sample_cmdline,
    project_df_file=project_df)

def cmdline():
    args = parsers['main'].parse_args()

    # get the function to be run
    if 'func' in args:
        func = args.func
    else:
        if args.main is not None:
            parsers[args.main].print_help()
        else:
            parsers['main'].print_help()
        return 0
        
    # get the args and delete the func key, get only set values
    args = {key: value for key, value in vars(args).items() if value is not None}
    args.pop('func', None)

    func(args)

if __name__ == "__main__":
    cmdline()
