import os
import yaml
import argparse
import re
import logging

from spacemake.errors import *
from spacemake.util import str2bool, assert_file, bool_in_str, message_aggregation
from spacemake.util import check_star_index_compatibility

logger_name = 'spacemake.config'

def get_puck_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False, add_help=False)
    parser.add_argument('--name',
        help='name of the puck', type=str, required=True)
    parser.add_argument('--width_um',
        type=float, required=required,
        help='width of the puck in microns')
    parser.add_argument('--spot_diameter_um',
        type=float, required=required,
        help='diameter of the spots in this puck, in microns')
    parser.add_argument('--barcodes',
        type=str, required=False,
        help='path to barcode file. of not provided the --puck_barcode_file variable' +
            ' of `spacemake projects add_sample` has to be set')

    return parser

def get_run_mode_parser(required=True):
    parser = argparse.ArgumentParser(
            allow_abbrev=False, formatter_class = argparse.RawTextHelpFormatter,
            description='add/update run_mode parent parser',
            add_help = False)
    parser.add_argument('--name', type=str,
            help='name of the run_mode to be added',
            required = True)
    parser.add_argument('--parent_run_mode', type=str,
            help='Name of the parent run_mode. All run_modes will fall back to \'default\'')
    parser.add_argument('--umi_cutoff', type=int, nargs='+',
            help='umi_cutoff for this run_mode.' +\
                    'the automated analysis will be run with these cutoffs')
    parser.add_argument('--n_beads',
            type=int,
            help='number of expected beads for this run_mode')
    parser.add_argument(
            '--clean_dge',
            required=False,
            choices=bool_in_str,
            type=str,
            help='if True, the DGE will be cleaned of barcodes which overlap with primers')
    parser.add_argument(
            '--detect_tissue',
            required=False,
            choices=bool_in_str,
            type=str,
            help='By default only beads having at least umi_cutoff UMI counts are analysed '+\
                    'during the automated analysis, all other beads are filtered out. If this ' +\
                    'parameter is set, contiguous islands within umi_cutoff passing beads will '+\
                    'also be included in the analysis')
    parser.add_argument(
            '--polyA_adapter_trimming',
            required=False,
            choices=bool_in_str,
            type=str,
            help='if set, reads will have polyA stretches and adapter sequence overlaps trimmed '+\
                    'BEFORE mapping.')
    parser.add_argument(
        '--count_intronic_reads',
        required=False,
        choices=bool_in_str,
        type=str,
        help='if set, INTRONIC reads will also be countsed (apart from UTR and CDS)')
    parser.add_argument('--count_mm_reads', required=False,
        choices=bool_in_str, type=str,
        help='if True, multi-mappers will also be counted. For every '+
            'multimapper only reads which have one unique read mapped' +
            'to a CDS or UTR region will be counted')
    
    parser.add_argument('--mesh_data', required=False,
        choices=bool_in_str, type=str,
        help ='if True, this data will be \'mehsed\': a hexagonal structured '+
            'meshgrid will be created, where each new spot will have diameter'+
            ' of --mesh_spot_diameter_um micron diameter and the spots will ' +
            'be spaced --mesh_spot_distance_um microns apart')
    parser.add_argument('--mesh_type', required=False,
        choices=['circle', 'hexagon'], type = str,
        help = 'circle: circles with diameter of --mesh_spot_diameter_um will be placed' +
                ' in a hex grid, where he distance between any two circles will be ' +
                '--mesh_spot_distance_um\n' +
                'hexagon: a mesh of touching hexagons will be created, with their ' +
                'centers being --mesh_spot_distance_um apart. This will cover all ' +
                'the data without any holes')
    parser.add_argument('--mesh_spot_diameter_um',
        type=float, required=False,
        help='diameter of mesh spot, in microns. to create a visium-style '+
            'mesh, use 55um')
    parser.add_argument('--mesh_spot_distance_um',
        type=float, required=False,
        help='distance between mesh spots in um. to create a visium-style '+
            'mesh use 100um')

    return parser

def get_species_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False, add_help=False)
    parser.add_argument('--name',
        help = 'name of the species to be added',
        type = str, required=True)
    parser.add_argument('--genome',
        help = 'path to the genome (.fa) file for the species to be added',
        type = str, required=required)
    parser.add_argument('--annotation',
        help = 'path to the annotation (.gtf) file for the species to be added',
        type = str, required=required)
    parser.add_argument('--rRNA_genome',
        help = 'path to the ribosomal-RNA genome (.fa) file for the species to be added',
        default=None,
        type = str)
    parser.add_argument('--STAR_index_dir',
        help = 'path to STAR index directory',
        type = str,
        required = False)

    return parser

def get_barcode_flavor_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='add/update barcode_flavor',
        add_help = False)
    parser.add_argument('--name',
        help = 'name of the barcode flavor', type = str,
        required = True)
    parser.add_argument('--umi',
        help = 'structure of UMI, using python\'s list syntax. Example: to set UMI to ' +
        '13-20 NT of Read1, use --umi r1[12:20]. It is also possible to use the first 8nt of '+
        'Read2 as UMI: --umi r2[0:8]',
        type=str, required=required)
    parser.add_argument('--cell_barcode',
        help = 'structure of CELL BARCODE, using python\'s list syntax. Example: to set'+
            ' the cell_barcode to 1-12 nt of Read1, use --cell_barcode r1[0:12]. It is also possible '+
            ' to reverse the CELL BARCODE, for instance with r1[0:12][::-1] (reversing the first 12nt of' +
            ' Read1, and assigning them as CELL BARCODE).',
        type=str, required=required)

    return parser

def get_variable_action_subparsers(config, parent_parser, variable):
    if variable == 'barcode_flavors':
        variable_singular = variable[:-1]
        variable_add_update_parser = get_barcode_flavor_parser
    elif variable == 'pucks':
        variable_singular = variable[:-1]
        variable_add_update_parser = get_puck_parser
    elif variable == 'run_modes':
        variable_singular = variable[:-1]
        variable_add_update_parser = get_run_mode_parser
    elif variable == 'species':
        variable_singular = variable
        variable_add_update_parser = get_species_parser
    

    command_help = {
        'list': f'list {variable} and their settings',
        'delete': f'delete {variable_singular}',
        'add': f'add a new {variable_singular}',
        'update': f'update an existing barcode_flavor'
    }

    # list command
    list_parser = parent_parser.add_parser(f'list_{variable}',
        description = command_help['list'],
        help = command_help['list'])
    list_parser.set_defaults(
        func=lambda args: list_variables_cmdline(config, args),
        variable=variable)

    func=lambda args: add_update_delete_variable_cmdline(config, args)

    # delete command
    delete_parser = parent_parser.add_parser(f'delete_{variable_singular}',
        description = command_help['delete'],
        help = command_help['delete'])
    delete_parser.add_argument('--name',
        help = f'name of the {variable_singular} to be deleted',
        type=str, required=True)
    delete_parser.set_defaults(
        func=func, action='delete',
            variable=variable)

    # add command
    add_parser = parent_parser.add_parser(f'add_{variable_singular}',
        parents=[variable_add_update_parser()],
        description = command_help['add'],
        help=command_help['add'])
    add_parser.set_defaults(
        func=func, action='add',
            variable=variable)

    # update command
    update_parser = parent_parser.add_parser(f'update_{variable_singular}',
        parents=[variable_add_update_parser(False)],
        description = command_help['update'],
        help = command_help['update'])
    update_parser.set_defaults(
        func=func, action='update',
            variable=variable)

def setup_config_parser(config, subparsers):
    parser_config = subparsers.add_parser('config', help = 'configure spacemake')
    parser_config_subparsers = parser_config.add_subparsers(help = 'config sub-command help')

    for variable in ConfigFile.main_variables_pl2sg.keys():
        get_variable_action_subparsers(
            config,
            parser_config_subparsers,
            variable)

    return parser_config

@message_aggregation(logger_name)
def add_update_delete_variable_cmdline(config, args):
    # set the name and delete from dictionary
    name = args['name']
    variable = args['variable']
    action = args['action']

    config.assert_main_variable(variable)

    # remove the args from the dict
    del args['action']
    del args['variable']
    del args['name']

    if action == 'add':
        func = config.add_variable
    elif action == 'update':
        func = config.update_variable
    elif action == 'delete':
        func = config.delete_variable

    var_variables = func(variable, name, **args)

    # print and dump config file
    config.logger.info(yaml.dump(var_variables, sort_keys=False))
    config.dump()

@message_aggregation(logger_name)
def list_variables_cmdline(config, args):
    variable = args['variable']
    del args['variable']

    config.logger.info(f'Listing {variable}')
    config.logger.info(yaml.dump(config.variables[variable]))

class ConfigMainVariable:
    def __init__(self, name, **kwargs):
        self.name = name
        self.variables = {}

        # loop over kwargs
        for variable_key, variable in kwargs.items():
            if variable_key not in self.variable_types:
                raise UnrecognisedConfigVariable(f'{variable_key}',
                    list(self.variable_types.keys()))
            elif self.variable_types[variable_key] == 'int_list':
                self.variables[variable_key] = [
                    int(x) for x in kwargs[variable_key]
                ]
            else:
                new_type = self.variable_types[variable_key]
                self.variables[variable_key] = new_type(kwargs[variable_key])

    def __str__(self):
        class_name = self.__class__.__name__
        return f'{class_name}: {self.name}. variables: {self.variables}'

    def update(self, other):
        self.variables.update(other.variables)


class RunMode(ConfigMainVariable):
    variable_types = {
        'parent_run_mode': str,
        'n_beads': int,
        'umi_cutoff': 'int_list',
        'clean_dge': bool,
        'detect_tissue': bool,
        'polyA_adapter_trimming': bool,
        'count_mm_reads': bool,
        'count_intronic_reads': bool,
        'mesh_data': bool,
        'mesh_type': str,
        'mesh_spot_diameter_um': int,
        'mesh_spot_distance_um': int
    }

    def has_parent(self):
        if 'parent_run_mode' in self.variables.keys():
            return True
        else:
            return False

    @property
    def parent_name(self):
        if self.has_parent():
            return self.variables['parent_run_mode']
        else:
            return None

class Puck(ConfigMainVariable):
    variable_types = {
        'barcodes': str,
        'spot_diameter_um': float,
        'width_um': int
    }

    @property
    def has_barcodes(self):
        return ('barcodes' in self.variables
                and self.variables['barcodes']
                and self.variables['barcodes'] != 'None')

class ConfigFile:
    initial_config_path = os.path.join(os.path.dirname(__file__),
        'data/config/config.yaml') 

    main_variables_pl2sg = {
        'pucks':'puck',
        'barcode_flavors':'barcode_flavor',
        'run_modes': 'run_mode',
        'species': 'species'
    }

    main_variables_sg2pl = {value: key for key, value in main_variables_pl2sg.items()}

    main_variable_sg2type = {
        'puck': str,
        'barcode_flavor': str,
        'run_mode': 'str_list',
        'species': str
    }

    def __init__(self):
        self.variables = {
            'root_dir': '',
            'temp_dir': '/tmp',
            'species': {},
            'barcode_flavors': {},
            'run_modes': {},
            'pucks': {},
        }
        self.file_path = 'config.yaml'
        self.logger = logging.getLogger(logger_name)

    @classmethod
    def from_yaml(cls, file_path):
        cf = cls()

        config_yaml_variables = None
        with open(file_path, 'r') as f:
            config_yaml_variables = yaml.load(f,
                    Loader=yaml.FullLoader)

        if config_yaml_variables is not None:
            cf.variables = config_yaml_variables

        cf.file_path = file_path

        if file_path != cf.initial_config_path:
            initial_config = ConfigFile.from_yaml(cf.initial_config_path)

            # correct variables to ensure backward compatibility
            cf.correct() 
            
            # check which variables do not exist, if they dont, 
            # copy them from initial config
            for main_variable in cf.main_variables_pl2sg:
                # update new main_variables
                if main_variable not in cf.variables:
                    cf.variables[main_variable] = initial_config.variables[main_variable]

            # deduce variables which have 'default' value. this is to ensure spacemake
            # always runs w/o errors downstream: ie when barcode flavor, run_mode or puck
            # is set to default
            cf.vars_with_default = [key for key, value in initial_config.variables.items()
                if 'default' in value]

            for var_with_default in cf.vars_with_default:
                default_val = initial_config.variables[var_with_default]['default']
                if 'default' not in cf.variables[var_with_default]:
                    cf.variables[var_with_default]['default'] = default_val
                else:
                    # update default run mode with missing values
                    default_val.update(cf.variables[var_with_default]['default'])
                    cf.variables[var_with_default]['default'] = default_val

        return cf

    def assert_main_variable(self, variable):
        if variable not in self.main_variables_pl2sg.keys():
            raise UnrecognisedConfigVariable(variable,
                list(self.main_variables_pl2sg.keys()))

    def correct(self):
        # ensures backward compatibility
        if 'pucks' not in self.variables:
            self.variables['pucks'] = {}
            if 'puck_data' in self.variables:
                if 'pucks' in self.variables['puck_data']:
                    self.variables['pucks'] = self.variables['puck_data']['pucks']
                    del self.variables['puck_data']['pucks']
        
        if 'barcode_flavors' not in self.variables:
            self.variables['barcode_flavors'] = {}
            if 'knowledge' in self.variables:
                if 'barcode_flavor' in self.variables['knowledge']:
                    self.variables['barcode_flavors'] = self.variables['knowledge']['barcode_flavor']
        
        if 'species' not in self.variables:
            # get all the species and save them in the right place
            # if species is empty, create a species dictionary
            self.variables['species'] = {}

            if 'knowledge' in self.variables:
                # extract all annotation info, if exists
                for species in self.variables['knowledge'].get('annotations', {}):
                    if species not in self.variables['species']:
                        self.variables['species'][species] = {}

                    self.variables['species'][species]['annotation'] = \
                        self.variables['knowledge']['annotations'][species]

                for species in self.variables['knowledge'].get('genomes', {}):
                    if species not in self.variables['species']:
                        self.variables['species'][species] = {}

                    self.variables['species'][species]['genome'] = \
                        self.variables['knowledge']['genomes'][species]

                for species in self.variables['knowledge'].get('rRNA_genomes', {}):
                    if species not in self.variables['species']:
                        self.variables['species'][species] = {}

                    self.variables['species'][species]['rRNA_genome'] = \
                        self.variables['knowledge']['rRNA_genomes'][species]

        if 'knowledge' in self.variables:
            del self.variables['knowledge']
        
        # correct run modes
        for run_mode_name, run_mode_variables in self.variables['run_modes'].items():
            variables = run_mode_variables.copy()
            for var in run_mode_variables:
                if not var in RunMode.variable_types:
                    del variables[var]

            self.variables['run_modes'][run_mode_name] = variables

    def dump(self):
        with open(self.file_path, 'w') as fo:
            fo.write(yaml.dump(self.variables))

    def set_file_path(self, file_path):
        self.file_path = file_path

    @property
    def puck_data(self):
        return self.variables['puck_data']

    def variable_exists(self, variable_name, variable_key):
        return variable_key in self.variables[variable_name]

    def assert_variable(self, variable_name, variable_key):
        self.assert_main_variable(variable_name)
        if not isinstance(variable_key, list):
            variable_key = [variable_key]

        for key in variable_key:
            if not self.variable_exists(variable_name, key):
                variable_name_sg = self.main_variables_pl2sg[variable_name]
                raise ConfigVariableNotFoundError(variable_name_sg,
                        key)

    def delete_variable(self, variable_name, variable_key):
        self.assert_variable(variable_name, variable_key)

        if variable_name in self.vars_with_default and variable_key == 'default':
            raise EmptyConfigVariableError(variable_name)

        variable_data = self.variables[variable_name][variable_key]

        del self.variables[variable_name][variable_key]

        return variable_data

    def process_run_mode_args(self, **kwargs):
        # typeset boolean values of run_mode
        default_run_mode = self.get_variable('run_modes', 'default')

        for key, value in default_run_mode.items():
            if isinstance(value, bool) and key in kwargs.keys():
                kwargs[key] = str2bool(kwargs[key])

        return kwargs

    def process_barcode_flavor_args(self, cell_barcode=None, umi=None):
        bam_tags = 'CR:{cell},CB:{cell},MI:{UMI},RG:{assigned}'

        # r(1|2) and then string slice
        to_match = r'r(1|2)(\[((?=-)-\d+|\d)*\:((?=-)-\d+|\d*)(\:((?=-)-\d+|\d*))*\])+$'

        if umi is not None and re.match(to_match, umi) is None:
            raise InvalidBarcodeStructure('umi', to_match)

        if cell_barcode is not None and re.match(to_match, cell_barcode) is None:
            raise InvalidBarcodeStructure('umi', to_match)

        barcode_flavor = {
                'bam_tags': bam_tags}

        if umi is not None:
            barcode_flavor['UMI'] = umi

        if cell_barcode is not None:
            barcode_flavor['cell'] = cell_barcode

        return barcode_flavor

    def process_species_args(
        self,
        genome=None,
        annotation=None,
        rRNA_genome = None,
        STAR_index_dir = None,
    ):
        assert_file(genome, default_value=None, extension = ['.fa', '.fa.gz'])
        assert_file(annotation, default_value=None, extension = ['.gtf', '.gtf.gz'])
        assert_file(rRNA_genome, default_value = None, extension = ['.fa'])

        species = {}

        if genome is not None:
            species['genome'] = genome
        
        if annotation is not None:
            species['annotation'] = annotation

        if rRNA_genome is not None:
            species['rRNA_genome'] = rRNA_genome

        if STAR_index_dir is not None:
            check_star_index_compatibility(STAR_index_dir)
            species['STAR_index_dir'] = STAR_index_dir

        return species

    def process_puck_args(self, width_um=None, spot_diameter_um=None, barcodes = None):
        assert_file(barcodes, default_value = None, extension = 'all')
        
        puck = {}

        if width_um is not None:
            puck['width_um'] = float(width_um),

        if spot_diameter_um is not None:
            puck['spot_diameter_um'] = float(spot_diameter_um)

        if barcodes is not None:
            puck['barcodes'] = barcodes

        return puck

    def process_variable_args(self, variable, **kwargs):
        if variable == 'barcode_flavors':
            return self.process_barcode_flavor_args(**kwargs)
        elif variable == 'run_modes':
            return self.process_run_mode_args(**kwargs)
        elif variable == 'pucks':
            return self.process_puck_args(**kwargs)
        elif variable == 'species':
            return self.process_species_args(**kwargs)
        else:
            ValueError(f'Invalid variable: {variable}')
        
    def add_variable(self, variable, name, **kwargs):
        if not self.variable_exists(variable, name):
            values = self.process_variable_args(variable, **kwargs)
            self.variables[variable][name] = values
        else:
            if variable in ['run_modes', 'pucks', 'barcode_flavors']:
                # drop the last s
                variable = variable[:-1]

            raise DuplicateConfigVariableError(variable, name)

        return kwargs

    def update_variable(self, variable, name, **kwargs):
        if self.variable_exists(variable, name):
            values = self.process_variable_args(variable, **kwargs)
            self.variables[variable][name].update(values)

            variable_data = self.variables[variable][name]
        else:
            if variable in ['run_modes', 'pucks', 'barcode_flavors']:
                # drop the last s
                variable = variable[:-1]

            raise ConfigVariableNotFoundError(variable, name)

        return variable_data

    def get_variable(self, variable, name):
        if not self.variable_exists(variable, name):
            raise ConfigVariableNotFoundError(variable, name)
        else:
            return self.variables[variable][name]

    def get_run_mode(self, name):
        # first load the default values
        rm = RunMode(name,
            **self.get_variable('run_modes', name))

        default_rm = RunMode('default',
            **self.get_variable('run_modes', 'default'))

        if rm.has_parent():
            parent_rm = self.get_run_mode(rm.parent_name)
            parent_rm.update(rm)
            rm = parent_rm

        # update with self
        default_rm.update(rm)
        rm = default_rm

        return rm

    def get_puck(self, name, return_empty=False):
        try:
            return Puck(name, **self.get_variable('pucks', name))
        except ConfigVariableNotFoundError as e:
            if not return_empty:
                raise
            else:
                return Puck(name)

