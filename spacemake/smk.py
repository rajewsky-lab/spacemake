import re
import sys
import os
import snakemake
import argparse
import yaml

from shutil import copyfile

config_path = 'config.yaml'

class ConfigFile:
    def __init__(self, file_path):
        self.file_path = file_path
        self.variables = yaml.load(open(file_path),
                    Loader=yaml.FullLoader)

    def dump(self):
        with open(self.file_path, 'w') as fo:
            fo.write(yaml.dump(self.variables))

    def set_file_path(self, file_path):
        self.file_path = file_path

    @classmethod
    def delete_run_mode(cls, args):
        name = args['name']
        cf = cls(config_path)
        if not name in cf.variables['run_modes'].keys():
            print(f'run_mode: {name} does not exists, so it cannot be deleted')
            return 0
        del cf.variables['run_modes'][name]
        cf.dump()

    @classmethod
    def list_run_modes(cls, args):
        cf = cls(config_path)

        run_modes = cf.variables['run_modes']

        print(yaml.dump(run_modes))

        return run_modes

    @classmethod
    def add_run_mode(cls, args):
        cf = cls(config_path)

        name = args['name']
        del args['name']
        
        # add and save new run mode
        if name in cf.variables['run_modes']:
            msg = f'run_mode: {name} already exists, so it cannot be added.\n'
            msg += 'you can update this run_mode using the `spacemake config update_run_mode`'
            msg += ' command.\nor delete it using the `spacemake config delete_run_mode`' 
            msg += ' command.'

            print(msg)
            return 0

        msg = f'adding run_mode: {name}\n'
        msg += '----------\n'
        msg += 'variables:\n'
        msg += '----------\n'
        msg += yaml.dump(args, sort_keys=False)

        print(msg)

        cf.variables['run_modes'][name] = args
        cf.dump()

        print(f'run_mode: {name} added')

        return 1
    
    @classmethod
    def update_run_mode(cls, args):
        cf = cls(config_path)

        name = args['name']
        del args['name']

        if not name in cf.variables['run_modes']:
            msg = f'run_mode: {name} does not exists, so it cannot be updated.\n'
            msg += 'add this run_mode using the `spacemake config add_run_mode` command'
            print(msg)
            return 0

        msg = f'updating run_mode: {name}\n'
        msg += '----------\n'
        msg += 'variables:\n'
        msg += '----------\n'
        msg += yaml.dump(args, sort_keys=False)

        print(msg)

        cf.variables['run_modes'][name].update(args)
        cf.dump()

        print(f'run_mode: {name} updated')

        return 1

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
parser_config = subparsers.add_parser('config', help = 'configure spacemake')
parser_config_subparsers = parser_config.add_subparsers(help = 'config sub-command help')


####################
# SPACEMAKE CONFIG #
####################
## list run_modes ##
parser_config_list_run_modes = parser_config_subparsers.add_parser('list_run_modes',
        help = 'list available run_modes')
parser_config_list_run_modes.set_defaults(func=ConfigFile.list_run_modes)

## add/update run_mode ##
# create the parent parser
parser_config_run_mode_parent = argparse.ArgumentParser(
    description='add/update run_mode parent parser',
    add_help = False)
parser_config_run_mode_parent.add_argument('--name', type=str, required=True,
    help='name of the run_mode to be added')
parser_config_run_mode_parent.add_argument('--parent_run_mode', type=str,
    help='Name of the parent run_mode. All run_modes will fall back to \'default\'')
parser_config_run_mode_parent.add_argument('--umi_cutoff', type=int, nargs='+',
    help='umi_cutoff for this run_mode.' +\
         'the automated analysis will be run with these cutoffs')
parser_config_run_mode_parent.add_argument(
    '--clean_dge',
    required=False,
    action=argparse.BooleanOptionalAction,
    help='if set, the DGE will be cleaned of barcodes which overlap with primers')
parser_config_run_mode_parent.add_argument(
    '--detect_tissue',
    required=False,
    action=argparse.BooleanOptionalAction,
    help='By default only beads having at least umi_cutoff UMI counts are analysed '+\
         'during the automated analysis, all other beads are filtered out. If this ' +\
         'parameter is set, contiguous islands within umi_cutoff passing beads will '+\
         'also be included in the analysis')
parser_config_run_mode_parent.add_argument(
    '--plot_bead_size',
    help='The bead size to be used when plotting the beads in 2D, during the '+\
         'automated report generation. Defaults to 1.')
parser_config_run_mode_parent.add_argument(
    '--polyA_adapter_trimming',
    required=False,
    action=argparse.BooleanOptionalAction,
    help='If set, reads will have polyA stretches and adapter sequence overlaps trimmed '+\
         'BEFORE mapping.')
parser_config_run_mode_parent.add_argument(
    '--count_intronic_reads',
    required=False,
    action=argparse.BooleanOptionalAction,
    help='If set, INTRONIC reads will also be countsed (apart from UTR and CDS)')
parser_config_run_mode_parent.add_argument(
    '--count_mm_reads',
    required=False,
    action=argparse.BooleanOptionalAction,
    help='If True, multi-mappers will also be counted. For every multimapper only reads which '+\
         'have one unique read mapped to a CDS or UTR region will be counted')
        
# add run_mode
parser_config_add_run_mode = parser_config_subparsers.add_parser('add_run_mode',
        help = 'add a new run_mode', parents=[parser_config_run_mode_parent])
parser_config_add_run_mode.set_defaults(func=ConfigFile.add_run_mode)

# update run mode
parser_config_update_run_mode = parser_config_subparsers.add_parser('update_run_mode',
    help = 'update run_mode', parents=[parser_config_run_mode_parent])
parser_config_update_run_mode.set_defaults(func=ConfigFile.update_run_mode)

parser_config_delete_run_mode = parser_config_subparsers.add_parser('delete_run_mode',
    help = 'delete a run_mode')
parser_config_delete_run_mode.add_argument('--name',
    required=True,
    type=str,
    help='run_mode to be deleted')
parser_config_delete_run_mode.set_defaults(func=ConfigFile.delete_run_mode)

def cmdline():
    args = parser.parse_args()

    # get the function to be run
    func = args.func
    
    # get the args and delete the func key, get only set values
    args = {key: value for key, value in vars(args).items() if value is not None}
    args.pop('func', None)

    func(args)

if __name__ == "__main__":
    cmdline()