import os
import yaml
import argparse
import re

from ..errors import *
from ._util import str2bool, assert_file, bool_in_str, LINE_SEPARATOR

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
        'plot_bead_size': float,
        'detect_tissue': bool,
        'polyA_adapter_trimming': bool,
        'count_mm_reads': bool,
        'count_intronic_reads': bool,
        'mesh_data': bool,
        'mesh_spot_diameter_um': int,
        'mesh_spot_distance_um': int
    }

    def has_parent(self):
        return 'parent_run_mode' in self.variables

    @property
    def parent_name(self):
        if self.has_parent():
            return self.variables['parent_run_mode']
        else:
            return None

class Puck(ConfigMainVariable):
    variable_types = {
        'barcodes': str,
        'spot_diameter_um': int,
        'width_um': int
    }

    @property
    def has_barcodes(self):
        return 'barcodes' in self.variables

class ConfigFile:
    initial_config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)),
        'config/config.yaml') 

    main_variables_pl2sg = {
        'pucks':'puck',
        'barcode_flavors':'barcode_flavor',
        'run_modes': 'run_mode',
        'species': 'species'
    }

    main_variables_sg2pl = {value: key for key, value in main_variables_pl2sg.items()}

    def __init__(self, file_path):
        self.file_path = file_path
        self.variables = yaml.load(open(file_path),
                    Loader=yaml.FullLoader)

        if file_path != self.initial_config_path:
            initial_config = ConfigFile.get_initial_config()

            # correct variables to ensure backward compatibility
            self.correct() 
            
            # check which variables do not exist, if they dont, 
            # copy them from initial config
            for main_variable in self.main_variables_pl2sg:
                # update new main_variables
                if main_variable not in self.variables:
                    self.variables[main_variable] = initial_config.variables[main_variable]

            # deduce variables which have 'default' value. this is to ensure spacemake
            # always runs w/o errors downstream: ie when barcode flavor, run_mode or puck
            # is set to default
            self.vars_with_default = [key for key, value in initial_config.variables.items()
                if 'default' in value]

            for var_with_default in self.vars_with_default:
                default_val = initial_config.variables[var_with_default]['default']
                if 'default' not in self.variables[var_with_default]:
                    self.variables[var_with_default]['default'] = default_val
                else:
                    # update default run mode with missing values
                    default_val.update(self.variables[var_with_default]['default'])
                    self.variables[var_with_default]['default'] = default_val

            if self.file_path != initial_config.file_path:
                # only dump if not initial config
                self.dump()

    @classmethod
    def get_initial_config(cls):
        cf = cls(cls.initial_config_path)

        return cf

    def assert_main_variable(self, variable):
        if variable not in self.main_variables_pl2sg.keys():
            raise UnrecognisedConfigVariable(variable,
                list(self.main_variables_pl2sg.keys()))

    def correct(self):
        # ensures backward compatibility
        if 'pucks' not in self.variables and 'pucks' in self.variables['puck_data']:
            self.variables['pucks'] = self.variables['puck_data']['pucks']
            del self.variables['puck_data']['pucks']

        if 'barcode_flavors' not in self.variables.keys() and 'barcode_flavor' in self.variables['knowledge']:
            self.variables['barcode_flavors'] = self.variables['knowledge']['barcode_flavor']
        
        if 'species' not in self.variables and 'knowledge' in self.variables:
            # get all the species and save them in the right place
            # if species is empty, create a species dictionary
            self.variables['species'] = {}

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
        
        if 'species' not in self.variables:
            self.variables['species'] = {}

    def dump(self):
        with open(self.file_path, 'w') as fo:
            fo.write(yaml.dump(self.variables))

    def set_file_path(self, file_path):
        self.file_path = file_path

    @property
    def puck_data(self):
        return self.variables['puck_data']

    def list_variables_cmdline(self, args):
        variable = args['variable']
        del args['variable']

        msg = f'Listing {variable}\n'
        msg += LINE_SEPARATOR
        msg += yaml.dump(self.variables[variable])

        print(msg)

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

    def __process_run_mode_args(self, **kwargs):
        # typeset boolean values of run_mode
        default_run_mode = self.get_variable('run_modes', 'default')

        for key, value in default_run_mode.items():
            if isinstance(value, bool) and key in kwargs.keys():
                kwargs[key] = str2bool(kwargs[key])

        return kwargs

    def __process_barcode_flavor_args(self, cell_barcode=None, umi=None):
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

    def __process_species_args(self, genome=None, annotation=None, rRNA_genome = None):
        assert_file(genome, default_value=None, extension = '.fa')
        assert_file(annotation, default_value=None, extension = '.gtf')
        assert_file(rRNA_genome, default_value = None, extension = '.fa')

        species = {}

        if genome is not None:
            species['genome'] = genome
        
        if annotation is not None:
            species['annotation'] = annotation

        if rRNA_genome is not None:
            species['rRNA_genome'] = rRNA_genome

        return species

    def __process_puck_args(self, width_um=None, spot_diameter_um=None, barcodes = None):
        assert_file(barcodes, default_value = None, extension = 'all')
        
        puck = {}

        if width_um is not None:
            puck['width_um'] = float(width_um),

        if spot_diameter_um is not None:
            puck['spot_diameter_um'] = float(spot_diameter_um)

        if barcodes is not None:
            puck['barcodes'] = barcodes

        return puck

    def __process_variable_args(self, variable, **kwargs):
        if variable == 'barcode_flavors':
            return self.__process_barcode_flavor_args(**kwargs)
        elif variable == 'run_modes':
            return self.__process_run_mode_args(**kwargs)
        elif variable == 'pucks':
            return self.__process_puck_args(**kwargs)
        elif variable == 'species':
            return self.__process_species_args(**kwargs)
        else:
            ValueError(f'Invalid variable: {variable}')
        
    def add_variable(self, variable, name, **kwargs):
        if not self.variable_exists(variable, name):
            values = self.__process_variable_args(variable, **kwargs)
            self.variables[variable][name] = kwargs
        else:
            if variable in ['run_modes', 'pucks', 'barcode_flavors']:
                # drop the last s
                variable = variable[:-1]

            raise DuplicateConfigVariableError(variable, name)

        return kwargs

    def update_variable(self, variable, name, **kwargs):
        if self.variable_exists(variable, name):
            values = self.__process_variable_args(variable, **kwargs)
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
        rm = RunMode('default',
            **self.get_variable('run_modes', 'default'))

        # update the default first
        rm.update(RunMode(name,
            **self.get_variable('run_modes', name)))

        if rm.has_parent():
            rm.update(self.get_run_mode(rm.parent_name))
        else:
            return rm

    def get_puck(self, name, return_empty=False):
        try:
            return Puck(name, **self.get_variable('pucks', name))
        except ConfigVariableNotFoundError as e:
            if not return_empty:
                raise
            else:
                return Puck(name)

    def add_update_delete_variable_cmdline(self, args):
        # set the name and delete from dictionary
        name = args['name']
        variable = args['variable']
        action = args['action']

        self.assert_main_variable(variable)

        # remove the args from the dict
        del args['action']
        del args['variable']
        del args['name']

        if action == 'add':
            msg = 'Adding'
            succ_msg = 'added'
            func = self.add_variable
        elif action == 'update':
            msg = 'Updating'
            succ_msg = 'after update'
            func = self.update_variable
        elif action == 'delete':
            msg = 'Deleting'
            succ_msg = 'deleted'
            func = self.delete_variable

        msg = f'{msg} {variable}: {name}\n'
        msg += LINE_SEPARATOR

        try:
            var_variables = func(variable=variable, name=name, **args)

            msg += f'SUCCESS: {variable}={name} {succ_msg} successfully.\n'
            msg += LINE_SEPARATOR
            msg += f'{succ_msg} {variable} variables:\n'
            msg += yaml.dump(var_variables, sort_keys=False)
            self.dump()
        except SpacemakeError as e:
            msg += str(e)
        finally:
            print(msg)

    def __get_puck_parser(self, required=True):
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

    def __get_species_parser(self, required=True):
        parser = argparse.ArgumentParser(
            allow_abbrev=False, add_help=False)
        parser.add_argument('--name',
            help = 'name of the species to be added',
            type = str, required=True)
        parser.add_argument('--genome',
            help = 'path to the genome (.fa) file for the species to be added',
            type = str, required=True)
        parser.add_argument('--annotation',
            help = 'path to the annotation (.gtf) file for the species to be added',
            type = str, required=True)
        parser.add_argument('--rRNA_genome',
            help = 'path to the ribosomal-RNA genome (.fa) file for the species to be added',
            default=None,
            type = str)

        return parser

    def __get_barcode_flavor_parser(self, required=True):
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

    def __get_run_mode_parser(self, required=True):
        parser = argparse.ArgumentParser(
                allow_abbrev=False,
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
                '--plot_bead_size',
                help='The bead size to be used when plotting the beads in 2D, during the '+\
                        'automated report generation. Defaults to 1.', type=float)
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
        parser.add_argument('--mesh_spot_diameter_um',
            type=float, required=False,
            help='diameter of mesh spot, in microns. to create a visium-style '+
                'mesh, use 55um')
        parser.add_argument('--mesh_spot_distance_um',
            type=float, required=False,
            help='distance between mesh spots in um. to create a visium-style '+
                'mesh use 100um')

        return parser

    def __add_variable_subparsers(self, parent_parser, variable):
        if variable == 'barcode_flavors':
            variable_singular = variable[:-1]
            variable_add_update_parser = self.__get_barcode_flavor_parser
        elif variable == 'pucks':
            variable_singular = variable[:-1]
            variable_add_update_parser = self.__get_puck_parser
        elif variable == 'run_modes':
            variable_singular = variable[:-1]
            variable_add_update_parser = self.__get_run_mode_parser
        elif variable == 'species':
            variable_singular = variable
            variable_add_update_parser = self.__get_species_parser
        

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
            func=self.list_variables_cmdline, variable=variable)

        # delete command
        delete_parser = parent_parser.add_parser(f'delete_{variable_singular}',
            description = command_help['delete'],
            help = command_help['delete'])
        delete_parser.add_argument('--name',
            help = f'name of the {variable_singular} to be deleted',
            type=str, required=True)
        delete_parser.set_defaults(
            func=self.add_update_delete_variable_cmdline, action='delete',
                variable=variable)

        # add command
        add_parser = parent_parser.add_parser(f'add_{variable_singular}',
            parents=[variable_add_update_parser()],
            description = command_help['add'],
            help=command_help['add'])
        add_parser.set_defaults(
            func=self.add_update_delete_variable_cmdline, action='add',
                variable=variable)

        # update command
        update_parser = parent_parser.add_parser(f'update_{variable_singular}',
            parents=[variable_add_update_parser(False)],
            description = command_help['update'],
            help = command_help['update'])
        update_parser.set_defaults(
            func=self.add_update_delete_variable_cmdline, action='update',
                variable=variable)

    def get_subparsers(self, subparsers):
        parser_config = subparsers.add_parser('config', help = 'configure spacemake')
        parser_config_subparsers = parser_config.add_subparsers(help = 'config sub-command help')

        # add run_mode parser
        self.__add_variable_subparsers(parser_config_subparsers,
            'run_modes')

        # add barcode_flavor parser
        self.__add_variable_subparsers(parser_config_subparsers,
            'barcode_flavors')

        # add pucks parser
        self.__add_variable_subparsers(parser_config_subparsers,
            'pucks')

        # add species parser
        self.__add_variable_subparsers(parser_config_subparsers,
            'species')

        return parser_config

