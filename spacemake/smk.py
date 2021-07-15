import re
import sys
import os
import snakemake
import argparse

from shutil import copyfile

config_path = 'config.yaml'

class ConfigFile:
    def __init__(self, file_path):
        variables = yaml.load(open(file_path),
                    Loader=yaml.FullLoader)

def spacemake_config(args):
    pass

def spacemake_init(args):
    if os.path.isfile(config_path):
        msg = "spacemake has already been initiated in this directory.\n"
        msg += "use other commands to run and analyse samples."
        print(f"{msg}")
        return 0

    initial_config = os.path.join(os.path.dirname(__file__),'config/config.yaml')

    copyfile(initial_config, config_path)

def spacemake_run(args):
    if not os.path.isfile(config_path):
        msg = "spacemake has not been initalised yet.\n"
        msg += "please run `spacemake init` to start a new project"
        return 0

    # get the snakefile
    snakefile = os.path.join(os.path.dirname(__file__), 'snakemake/main.smk')
    # run snakemake
    snakemake.snakemake(snakefile, configfiles=[config_path],
            cores = args.cores, dryrun=args.dryrun,
            config={'root_dir': '', 'temp_dir': '/tmp'})

## define parser
parser = argparse.ArgumentParser(description='spacemake: bioinformatic pipeline for processing and analysis of spatial-transcriptomics data')

subparsers = parser.add_subparsers(help='sub-command help')

## spacemake_init args
parser_init = subparsers.add_parser('init', help = 'initialise spacemake')
parser_init.set_defaults(func=spacemake_init)

## spacemake_run args
parser_run = subparsers.add_parser('run', help = 'run spacemake')
parser_run.add_argument('--cores', type=int, default=1, help = 'number of cores to be used in total')
parser_run.add_argument('--dryrun', '-n', action='store_true', help = 'invokes a dry snakemake run, printing only commands')
parser_run.set_defaults(func=spacemake_run)

## spacemake_config args

def cmdline():
    args = parser.parse_args()

    args.func(args)

if __name__ == "__main__":
    cmdline()
