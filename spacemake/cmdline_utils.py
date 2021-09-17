import pandas as pd
import os
import errno
import yaml
import math
import argparse
import datetime
import re
from functools import reduce
from operator import getitem, setitem
from spacemake.errors import *

LINE_SEPARATOR = '-'*50+'\n'

bool_in_str = ['True', 'true', 'False', 'false']

def assert_file(file_path, default_value='none', extension='all'):
    if file_path == default_value:
        # file doesn't exist but has the default value,
        # so we do not need to assert anything
        return 0

    # check if file exists, raise error if not
    if not os.path.isfile(file_path):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), file_path)


    if not file_path.endswith(extension) and extension != 'all':
        raise FileWrongExtensionError(file_path, extension)

    return 0

def str2bool(var):
    if isinstance(var, bool):
        return var

    if var in ['True', 'true']:
        return True
    elif var in ['False', 'false']:
        return False
    else:
        raise ValueError(f'variable should be boolean, or one of: {bool_in_str}')
        
def set_dict_variable(dictionary, path, new_val):
    if len(path) == 1:
        dictionary[path[0]] = new_val
        return dictionary
    elif len(path) > 1:
        dictionary[path[0]] = set_dict_variable(dictionary[path[0]],
            path[1:], new_val)

        return dictionary
    else:
        raise ValueError(f'path cant be an empty list')

class ConfigFile:
    initial_config_path = os.path.join(os.path.dirname(__file__),
        'config/config.yaml') 

    main_variables_pl2sg = {
        'pucks':'puck',
        'barcode_flavors':'barcode_flavor',
        'run_modes': 'run_mode',
        'species': 'species'
    }

    main_variables_sg2pl = {value: key for main_variables_pl2sg.items()}

    def __init__(self, file_path):
        self.file_path = file_path
        self.variables = yaml.load(open(file_path),
                    Loader=yaml.FullLoader)

        if file_path != self.initial_config_path:
            initial_config = ConfigFile.get_initial_config()
            default_run_mode = initial_config.variables['run_modes']['default']

            if 'default' not in self.variables['run_modes'].keys():
                self.variables['run_modes']['default'] = default_run_mode
            else:
                # update default run mode with missing values
                default_run_mode.update(self.variables['run_modes']['default'])
                self.variables['run_modes']['default'] = default_run_mode

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
        if 'pucks' not in self.variables.keys():
            self.variables['pucks'] = self.variables['puck_data']['pucks']
            del self.variables['puck_data']['pucks']

        if 'barcode_flavors' not in self.variables.keys():
            self.variables['barcode_flavors'] = self.variables['knowledge']['barcode_flavor']
        
        if 'species' not in self.variables.keys():
            # get all the species and save them in the right place
            self.variables['species'] = {}

            # extract all annotation info
            for species in self.variables['knowledge'].get('annotations', {}).keys():
                if species not in self.variables['species'].keys():
                    self.variables['species'][species] = {}

                self.variables['species'][species]['annotation'] = \
                    self.variables['knowledge']['annotations'][species]

            for species in self.variables['knowledge'].get('genomes', {}).keys():
                if species not in self.variables['species'].keys():
                    self.variables['species'][species] = {}

                self.variables['species'][species]['genome'] = \
                    self.variables['knowledge']['genomes'][species]

            for species in self.variables['knowledge'].get('rRNA_genomes', {}).keys():
                if species not in self.variables['species'].keys():
                    self.variables['species'][species] = {}

                self.variables['species'][species]['rRNA_genome'] = \
                    self.variables['knowledge']['rRNA_genomes'][species]

        if 'knowledge' in self.variables.keys():
            del self.variables['knowledge']


        self.dump()


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

    def variable_exists(self, variable, name):
        return name in self.variables[variable].keys()

    def assert_variable(self, variable, name):
        if not self.variable_exists(variable, name):
            raise ConfigVariableNotFoundError(variable, name)

    def delete_variable(self, variable, name):
        if not self.variable_exists(variable, name):
            if variable in ['run_modes', 'pucks', 'barcode_flavors']:
                # drop the last s
                variable = variable[:-1]

            raise ConfigVariableNotFoundError(variable, name)
        else:
            variable_data = self.variables[variable][name]

            del self.variables[variable][name]

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

class ProjectDF:
    # default values of the project dataframe columns
    project_df_default_values = {
        "puck_id": "no_optical_puck",
        "sample_sheet": "none",
        "species": "none",
        "demux_barcode_mismatch": 1,
        "demux_dir": "none",
        "basecalls_dir": "none",
        "R1": "none",
        "R2": "none",
        "investigator": "unknown",
        "sequencing_date": "unknown",
        "experiment": "unknown",
        "puck_barcode_file": "none",
        "run_mode": ["default"],
        "barcode_flavor": "default",
        "is_merged":False,
        "merged_from":[],
        "puck":"default"}
    
    def __init__(
        self,
        file_path,
        config: ConfigFile = None
    ):
        self.file_path = file_path
        self.config = config

        if os.path.isfile(file_path):
            df = pd.read_csv(file_path,
                index_col=['project_id', 'sample_id'],
                converters={'run_mode': eval, 'merged_from': eval})
            project_list = []

            # update with new columns, if they exist.
            for ix, row in df.iterrows():
                s = pd.Series(self.project_df_default_values)
                s.update(row)
                s.name = row.name
                project_list.append(s)

            self.df = pd.concat(project_list, axis=1).T
            self.df.index.names = ['project_id', 'sample_id']

            # dump the result
            self.dump()
        else:
            index = pd.MultiIndex(
                names=['project_id', 'sample_id'],
                levels=[[],[]],
                codes=[[],[]]) 
            self.df = pd.DataFrame(columns = self.project_df_default_values.keys(),
                index=index)

    def __compute_max_barcode_mismatch(self, indices):
        """computes the maximum number of mismatches allowed for demultiplexing based
        on the indices present in the sample sheet."""
        num_samples = len(indices)

        if num_samples == 1:
            return 4
        else:
            max_mismatch = 3
            for i in range(num_samples - 1):
                for j in range(i + 1, num_samples):
                    hd = self.__hamming_distance(indices[i], indices[j])
                    max_mismatch = min(max_mismatch, math.ceil(hd / 2) - 1)
        return max_mismatch

    def __hamming_distance(self, string1, string2):
        return sum(c1 != c2 for c1, c2 in zip(string1, string2))
    
    def __find_barcode_file(self, puck_id):
        # first find directory of puck file

        # return none or the path of the file
        def get_barcode_file(path):
            if os.path.isfile(path):
                return path

            return "none"

        def find_dir(name, path):
            for root, dirs, files in os.walk(path):
                if name in dirs:
                    return os.path.join(root, name)

        puck_dir = find_dir(puck_id, self.config.puck_data['root'])
        path = None

        if puck_dir is not None:
            # puck dir exists, look for barcode file pattern
            path = os.path.join(puck_dir, self.config.puck_data["barcode_file"])

            return get_barcode_file(path)
        else:
            return self.project_df_default_values['puck_barcode_file']


    def get_sample_info(self, project_id, sample_id):
        # returns sample info from the projects df
        out_dict = self.df.loc[(project_id, sample_id)].to_dict()

        return out_dict

    def is_spatial(self, sample_id, project_id):
        puck_barcode_file = self.get_metadata('puck_barcode_file',
            sample_id = sample_id,
            project_id = project_id)

        puck = self.get_metadata('puck',
            project_id = project_id,
            sample_id = sample_id)

        puck_has_barcodes = False

        if puck != self.project_df_default_values['puck']:
            if 'barcodes' in self.config.get_puck(puck):
                puck_type_has_barcodes = True
        
        if puck_barcode_file != 'none' or puck_type_has_barcodes:
            return True
        else:
            return False

    def get_puck(self, project_id, sample_id):
        puck = self.get_metadata('puck',
            project_id = project_id,
            sample_id = sample_id)

        return self.config.get_variable('pucks', name=puck)

    def get_metadata(self, field, sample_id=None, project_id=None, **kwargs):
        df = self.df
        if sample_id is not None:
            df = df.query('sample_id == @sample_id')

        if project_id is not None:
            df = df.query('project_id == @project_id')

        for key, value in kwargs.items():
            df = df.loc[df.loc[:, key] == value]

        return df[field].to_list()[0]

    def dump(self):
        self.df.to_csv(self.file_path)

    def add_sample_sheet(self, sample_sheet_path, basecalls_dir):
        with open(sample_sheet_path) as sample_sheet:
            ix = 0
            investigator = "none"
            sequencing_date = "none"

            for line in sample_sheet:
                line = line.strip("\n")
                if "Investigator" in line:
                    investigator = line.split(",")[1]
                if "Date" in line:
                    sequencing_date = line.split(",")[1]
                if "[Data]" in line:
                    # the counter ix stops here
                    break
                else:
                    ix = ix + 1

        # read everything after [Data]
        df = pd.read_csv(sample_sheet_path, skiprows=ix + 1)
        # rename columns
        to_rename={
                "Sample_ID": "sample_id",
                "Sample_Name": "puck_id",
                "Sample_Project": "project_id",
                "Description": "experiment",
                "index": "index"
            }
        df.rename(
            columns=to_rename,
            inplace=True,
        )
        # select only renamed columns
        df = df[to_rename.values()]
        df["species"] = df["experiment"].str.split("_").str[-1]
        df["investigator"] = investigator
        df["sequencing_date"] = sequencing_date

        # rename columns 
        df["basecalls_dir"] = basecalls_dir
        df["demux_barcode_mismatch"] = self.__compute_max_barcode_mismatch(df["index"])
        df["sample_sheet"] = sample_sheet_path
        df["demux_dir"] = df["sample_sheet"].str.split("/").str[-1].str.split(".").str[0]
        df["puck_barcode_file"] = df.puck_id.apply(self.__find_barcode_file)
        df.set_index(['project_id', 'sample_id'], inplace=True)

        for ix, row in df.iterrows():
            self.__add_update_sample(ix[0], ix[1], **row.to_dict())

    def sample_exists(self, project_id = None, sample_id = None):
        if project_id is None or sample_id is None:
            raise Exception(f'you need to provide a sample_id and project_id to check if sample exists')
        else:
            ix = (project_id, sample_id) 

            return ix in self.df.index

    def __add_update_sample(self, project_id = None,
            sample_id = None, return_series = False,
            **kwargs):
        """
        adds or updates a sample with a given project_id and sample_id
        """
        ix = (project_id, sample_id) 

        # assert files first
        assert_file(
            kwargs.get('R1',
                       self.project_df_default_values['R1']),
            default_value=self.project_df_default_values['R1'],
            extension = '.fastq.gz')

        assert_file(
            kwargs.get('R2',
                       self.project_df_default_values['R2']),
            default_value=self.project_df_default_values['R2'],
            extension = '.fastq.gz')

        assert_file(
            kwargs.get('puck_barcode_file',
                        self.project_df_default_values['puck_barcode_file']),
            default_value=self.project_df_default_values['puck_barcode_file'])

        # check if run mode exists
        for run_mode in kwargs.get('run_mode', []):
            if not self.config.variable_exists('run_modes', run_mode):
                raise ConfigVariableNotFoundError('run_modes', run_mode)

        config_variables_to_check = {
            'pucks': 'puck',
            'barcode_flavors': 'barcode_flavor',
            'species':'species'}

        for cv_plural, cv_singular in config_variables_to_check.items():
            if cv_singular in kwargs.keys():
                if not self.config.variable_exists(cv_plural, kwargs[cv_singular]):
                    raise ConfigVariableNotFoundError(cv_singular,
                            kwargs[cv_singular])

        if self.sample_exists(project_id, sample_id):
            new_project = self.df.loc[ix].copy()
            new_project.update(pd.Series(kwargs))
            self.df.loc[ix] = new_project
        else:
            new_project = pd.Series(self.project_df_default_values)
            new_project.name = ix
            new_project.update(kwargs)

            self.df = self.df.append(new_project)

        if return_series:
            return new_project

    def add_sample(self, project_id, sample_id, **kwargs):
        if self.sample_exists(project_id, sample_id):
            raise SampleAlreadyExistsError((project_id, sample_id))
        else:
            return self.__add_update_sample(project_id, sample_id,
                return_series=True, **kwargs)

    def update_sample(self, project_id, sample_id, **kwargs):
        if self.sample_exists(project_id, sample_id):
            return self.__add_update_sample(project_id, sample_id,
                return_series=True, **kwargs)
        else:
            raise ProjectSampleNotFoundError('(project_id, sample_id)',
                (project_id, sample_id))

    def delete_sample(self, project_id, sample_id):
        if self.sample_exists(project_id, sample_id):
            ix = (project_id, sample_id)
            element = self.df.loc[ix]
            self.df.drop(ix, inplace=True)
            return element
        else:
            raise ProjectSampleNotFoundError('(project_id, sample_id)',
                (project_id, sample_id))

    def add_samples_from_yaml(self, projects_yaml_file):
        config = yaml.load(open(projects_yaml_file),
                Loader=yaml.FullLoader)
        demux_projects = config.get('projects', None)
        barcode_flavors = config.get('barcode_flavor', None)

        if demux_projects is not None:
            # if we have projects in the config file
            # get the samples
            for ip in demux_projects:
                self.add_sample_sheet(ip['sample_sheet'], ip['basecalls_dir'])

        # add additional samples from config.yaml, which have already been demultiplexed.
        for project in config['additional_projects']:
            self.__add_update_sample(**project)

        for barcode_flavor, to_set in barcode_flavors.items():
            projects = to_set.get('projects', [])
            samples = to_set.get('samples', [])

            print(projects, samples)
            self.set_barcode_flavor(barcode_flavor, projects, samples)


    def __get_project_sample_parser(self):
        parser = argparse.ArgumentParser(
            allow_abbrev=False,
            add_help=False)

        parser.add_argument('--project_id',
            type=str, required=True,
            help = 'project_id of the sample to be added/updated')

        parser.add_argument('--sample_id',
            type=str, required=True,
            help = 'sample_id of the sample to be added/updated')

        return parser

    def __get_sample_extra_arguments_parser(self):
        parser = argparse.ArgumentParser(
            allow_abbrev=False,
            add_help=False)

        parser.add_argument('--puck_id',
            type=str,
            help = 'puck_id of the sample to be added/update')

        parser.add_argument('--puck_barcode_file',
            type=str,
            help = 'the path to the file contining (x,y) positions of the barcodes')

        parser.add_argument('--investigator',
            type=str,
            help = 'add the name of the investigator(s) responsible for this sample')

        parser.add_argument('--experiment',
            type=str,
            help='description of the experiment')

        parser.add_argument('--sequencing_date',
            type=datetime.date.fromisoformat,
            help = 'sequencing date of the sample. format: YYYY-MM-DD')

        parser.add_argument('--run_mode',
            type=str,
            nargs='+',
            help = 'run_mode names for this sample.\n' +\
                'the sample will be processed using the provided run_modes.\n' +\
                'for merged samples, if left empty, the run_modes of the \n' +\
                'merged (input) samples will be intersected.\n')
        parser.add_argument('--puck', type=str,
            help = 'name of the puck for this sample. if puck contains a '+
                '`barcodes` path to a coordinate file, those coordinates' +
                ' will be used when processing this sample. if ' +
                ' not provided, a default puck will be used with '+
                'width_um=3000, spot_diameter_um=10')

        return parser

    def __get_barcode_flavor_species_parser(self, species_required=False):
        parser = argparse.ArgumentParser(
            allow_abbrev=False,
            add_help = False)

        parser.add_argument('--barcode_flavor',
            type=str,
            help = 'barcode flavor for this sample')

        parser.add_argument('--species',
            type=str,
            help = 'add the name of the species for this sample',
            required=species_required)

        return parser

    def __get_read_parser(self, reads_required=False):
        parser = argparse.ArgumentParser(
            allow_abbrev=False,
            add_help = False)

        parser.add_argument('--R1',
            type=str,
            help = '.fastq.gz file path to R1 reads',
            required=reads_required)

        parser.add_argument('--R2',
            type=str,
            help = '.fastq.gz file path to R2 reads',
            required=reads_required)

        return parser
    
    def add_update_delete_sample_cmdline(self, args):
        project_id = args['project_id']
        sample_id = args['sample_id']
        action = args['action']

        # remove the action from args
        del args['action']

        if action == 'add':
            msg = 'Adding'
            succ_msg = 'added'
            func = self.add_sample
        elif action == 'update':
            msg = 'Updating'
            succ_msg = 'updated'
            func = self.update_sample
        elif action == 'delete':
            msg = 'Deleting'
            succ_msg = 'deleted'
            func = self.delete_sample

        msg = f'{msg} sample: ({project_id}, {sample_id})\n'
        msg += LINE_SEPARATOR

        try:
            sample = func(**args)

            msg += f'SUCCESS: sample {succ_msg} successfully.\n'
            msg += LINE_SEPARATOR
            msg += str(sample)
            self.dump()
        except (FileNotFoundError, SpacemakeError) as e:
            msg += str(e)
        finally:
            print(msg)

    def __get_add_sample_sheet_parser(self):
        parser = argparse.ArgumentParser(
            allow_abbrev=False,
            description = 'add a new sample sheet to the samples',
            add_help=False)

        parser.add_argument('--sample_sheet', type = str,
            help = 'the path to the Illumina sample sheet',
            required=True)
        parser.add_argument('--basecalls_dir', type = str,
            help = 'path to the basecalls directory',
            required=True)

        return parser

    def add_sample_sheet_cmdline(self, args):
        self.add_sample_sheet(args['sample_sheet'],
            args['basecalls_dir'])

        self.dump()

    def add_samples_from_yaml_cmdline(self, args):
        self.add_samples_from_yaml(args['samples_yaml'])

        self.dump()

    def __get_project_sample_lists_parser(self, help_extra=''):
        parser = argparse.ArgumentParser(add_help=False)
        allow_abbrev=False,
        parser.add_argument('--project_id_list',
            nargs='*',
            default=[],
            type=str,
            help = 'project id-s for which we should ' + help_extra)
        parser.add_argument('--sample_id_list',
            nargs='*',
            default=[],
            type=str,
            help = 'sample id-s for which we should ' + help_extra)

        return parser

    def __set_species(self, species, project_id_list, sample_id_list):
        # raise error if both lists are empty
        if project_id_list == [] and sample_id_list == []:
            raise NoProjectSampleProvidedError()

        self.__assert_projects_samples_exist(project_id_list, sample_id_list)

        ix = self.df.query(
            'project_id in @projects or sample_id in @samples').index

        if self.config.variable_exists('species', species):
            self.df.loc[ix, 'species'] = species

            return ix.to_list()
        else:
            raise ConfigVariableNotFoundError('species', species)

    def set_species_cmdline(self, args):
        projects = args['project_id_list']
        samples = args['sample_id_list']
        species = args['species_name']

        msg = ''

        try:
            set_indices = self.__set_species(species, projects, samples)

            msg += f'Setting species: {species} for the'
            msg += f' following (project_id, sample_id) samples:\n'
            msg += f'{set_indices}\n'
            msg += LINE_SEPARATOR
            msg += 'SUCCESS: species set successfully.\n'

            self.dump()
        except (NoProjectSampleProvidedError, ConfigVariableNotFoundError,
                ProjectSampleNotFoundError) as e:
            msg += str(e)
        finally:
            print(msg)


    def __set_variable(self, ix, variable_name, variable_value, keep_old=False):
        variable_name_pl = self.config.main_variable_sg2pl[variable_name]
        self.config.assert_main_variable(variable_name_pl)
        self.config.assert_variable(variable_name_pl, variable_value)
        
        # add to current variable
        i_variable = self.df.at[ix, 'run_mode']
        i_run_mode.append(run_mode)
        self.df.at[ix, 'run_mode'] = list(set(i_run_mode))

    def __remove_variable(self, ix, variable_name, variable_value):
        variable_name_pl = self.config.main_variable_sg2pl[variable_name]
        self.config.assert_main_variable(variable_name_pl)
        self.config.assert_variable(variable_name_pl, variable_value)

        new_variable = [rm for rm in self.df.at[ix, 'run_mode'] if rm != run_mode]
        self.df.at[ix, 'run_mode'] = i_run_mode

    def __set_remove_variable(self, variable_name, variable_value,
            action, project_id_list = [], sample_id_list = [], keep_old = False):
        # raise error if both lists are empty
        if project_id_list == [] and sample_id_list == []:
            raise NoProjectSampleProvidedError()

        self.__assert_projects_samples_exist(projects, samples)

        variable_name_pl = self.config.main_variable_sg2pl[variable_name]
        self.config.assert_main_variable(variable_name_pl)

        # of only one provided use that, if both use intersection
        if project_id_list == []:
            ix = self.df.query(
                'sample_id in @sample_id_list').index
        elif sample_id_list == []:
            ix = self.df.query(
                'sample_id in @project_id_list').index
        else:
            ix = self.df.query(
                'project_id in @project_id_list and sample_id in @sample_id_list').index
            'project_id in @projects or sample_id in @samples').index

        if isinstance(variable_value, list):
            for var in variable_value:
                self.config.assert_variable(variable_name_pl, var)
        else:
            self.config.assert_variable(variable_name_pl, variable_value)
        
        for i, row in self.df.loc[ix, :].iterrows():
            # add/remove variable. if it already exists dont do anything
            if action == 'set':
                self.__set_variable(i, variable_name, variable_value, keep_old)
            elif action == 'remove':
                self.__remove_run_mode(i, run_mode)

        return ix.to_list()

    def set_remove_run_mode_cmdline(self, args):
        projects = args['project_id_list']
        samples = args['sample_id_list']
        run_modes = args['run_mode']
        action = args['action']

        if action == 'set':
            succ_msg = 'set'
            action_msg = 'Setting'
        elif action == 'remove':
            succ_msg = 'removed'
            action_msg = 'Removing'

        msg =''

        try:
            set_indices = self.__set_remove_run_mode(run_modes, action, projects, samples)
             
            msg += f'{action_msg} {run_modes} for the following (project_id, sample_id)'
            msg += f' samples:\n{set_indices}\n'
            msg += LINE_SEPARATOR
            msg += f'SUCCESS: run mode: {run_modes} {succ_msg} succesfully.\n'

            self.dump()
        except SpacemakeError as e:

            msg += str(e)
        finally:
            print(msg)

    def list_projects_cmdline(self, args):
        projects = args['project_id_list']
        samples = args['sample_id_list']
        variables = args['always_show'] + args['variables']

        df = self.df

        msg = 'Listing project/sample info\n'

        if projects != [] or samples != []:
            df = df.query('project_id in @projects or sample_id in @samples')
            msg += f'for projects: {projects} and samples: {samples}\n'
        else:
            msg += 'for all projects and samples\n'

        msg += LINE_SEPARATOR
        msg += f'Variables used: {variables}\n'
        msg += LINE_SEPARATOR
        
        # print the table
        msg += df.loc[:, variables].__str__()

        print(msg)

    def __assert_projects_samples_exist(self,
            project_id_list = [], sample_id_list = []):
        for project in project_id_list:
            if project not in self.df.index.get_level_values('project_id'):
                raise ProjectSampleNotFoundError('project_id', project)

        for sample in sample_id_list:
            if sample not in self.df.index.get_level_values('sample_id'):
                raise ProjectSampleNotFoundError('sample_id', sample)
        

    def __set_barcode_flavor(self, barcode_flavor,
            project_id_list = [], sample_id_list = []):
        # raise error if both lists are empty
        if project_id_list == [] and sample_id_list == []:
            raise NoProjectSampleProvidedError()

        self.__assert_projects_samples_exist(project_id_list, sample_id_list)

        if self.config.variable_exists('barcode_flavors', barcode_flavor):
            ix = self.df.query(
                'project_id in @project_id_list or sample_id in @sample_id_list').index

            if ix is None:
                raise NoProjectSampleProvidedError()

            self.df.loc[ix, 'barcode_flavor'] = barcode_flavor

            return ix.to_list()
        else:
            raise ConfigVariableNotFoundError('barcode_flavor', barcode_flavor)

    def set_barcode_flavor_cmdline(self, args):
        barcode_flavor = args['barcode_flavor']
        
        msg = ''

        try:
            set_indices = self.__set_barcode_flavor(**args)

            msg += f'Setting barcode flavor: {barcode_flavor} for the'
            msg += f' following (project_id, sample_id) samples:\n'
            msg += f'{set_indices}\n'
            msg += LINE_SEPARATOR
            msg += 'SUCCESS: barcode_flavor set successfully.\n'
            
            self.dump()

        except (ConfigVariableNotFoundError, NoProjectSampleProvidedError,
            ProjectSampleNotFoundError) as e:
            msg += str(e)
        finally:
            print(msg)

    def merge_samples(self,
        merged_project_id,
        merged_sample_id,
        project_id_list = [],
        sample_id_list = [],
        **kwargs
    ):
        if project_id_list == [] and sample_id_list == []:
            raise NoProjectSampleProvidedError()

        # check if projects and samples with these IDs exist
        self.__assert_projects_samples_exist(project_id_list, sample_id_list)

        if project_id_list == []:
            ix = self.df.query(
                'sample_id in @sample_id_list').index
        elif sample_id_list == []:
            ix = self.df.query(
                'sample_id in @project_id_list').index
        else:
            ix = self.df.query(
                'project_id in @project_id_list and sample_id in @sample_id_list').index

        consistent_variables = ['species', 'barcode_flavor', 'puck']

        ix_list = ix.to_list()

        if ix_list == []:
            raise ProjectSampleNotFoundError(
                '(project_id_list, sample_id_list)',
                (project_id_list, sample_id_list))

        # check for variable inconsistency
        # raise error if variable different between samples
        for variable in consistent_variables:
            variable_val = self.df.loc[ix, variable].to_list()

            if len(set(variable_val)) > 1:
                raise InconsistentVariablesDuringMerge(variable, ix_list)
            else:
                # attach the deduced, consisten variable
                kwargs[variable] = variable_val[0]

        
        variables_to_deduce = [
            'investigator',
            'experiment',
            'sequencing_date'
        ]

        for variable in variables_to_deduce:
            if variable not in kwargs.keys():
                kwargs[variable] = ';'.join(self.df.loc[ix, variable].unique())

            
        # if no run_mode provided, overwrite with user defined one
        if 'run_mode' not in kwargs.keys():
            run_mode = [set(run_mode)
                for run_mode in self.df.loc[ix].run_mode.to_list()
            ]

            # join run modes from parent samples
            if len(run_mode) == 1:
                run_mode = run_mode[0]
            else:
                run_mode = run_mode[0].intersection(*run_mode[1:])

            # create a list from the set intersection
            run_mode = list(run_mode)

            # if there are no common elements, throw an error
            if len(run_mode) == 0:
                raise InconsistentVariablesDuringMerge('run_mode', ix_list)
            
            # finally add run mode to arguments
            kwargs['run_mode'] = run_mode
        
        sample_added = self.add_sample(
            project_id = merged_project_id,
            sample_id = merged_sample_id,
            is_merged=True,
            merged_from = ix.to_list(),
            **kwargs)

        return (sample_added, ix)

    def merge_samples_cmdline(self, args):
        msg = ''
        try:
            sample, set_indices = self.merge_samples(**args)

            set_indices = set_indices.to_list()

            msg += f'Merging samples {set_indices}.\n'
            msg += LINE_SEPARATOR
            msg += 'SUCCESS: samples merged successfully.\n'
            msg += LINE_SEPARATOR
            msg += str(sample)

            self.dump()
        except (NoProjectSampleProvidedError, ConfigVariableNotFoundError,
                ProjectSampleNotFoundError, InconsistentVariablesDuringMerge) as e:
            msg += str(e)
        finally:
            print(msg)


    def get_subparsers(self, subparsers):
        parser = subparsers.add_parser('projects',
                help='manage projects and samples',
                description ='Using one of the subcommands specified below, it is possible to' +\
                        ' add/update/remove projects and their settings') 
        projects_subparser = parser.add_subparsers()

        help_desc = {
            'add_sample_sheet' : 'add projects and samples from Illumina sample sheet',
            'add_samples_from_yaml' : 'add several samples at once from a .yaml file',
            'add_sample': 'add a single sample from the command line',
            'update_sample':'update_existing_sample',
            'delete_sample':'delete existing sample',
            'merge_samples': 'merge several samples into one. ' +\
                    'samples need to have the same species',
            'list':'list several project(s) and sample(s) and their settings'
        }


        # ADD SAMPLE SHEET
        sample_add_sample_sheet = projects_subparser.add_parser('add_sample_sheet',
            description = help_desc['add_sample_sheet'],
            help = help_desc['add_sample_sheet'],
            parents=[self.__get_add_sample_sheet_parser()])
        sample_add_sample_sheet.set_defaults(func=self.add_sample_sheet_cmdline)

        # ADD SAMPLES FROM YAML
        sample_add_samples_yaml = projects_subparser.add_parser('add_samples_from_yaml',
            description = help_desc['add_samples_from_yaml'],
            help = help_desc['add_samples_from_yaml'])
        sample_add_samples_yaml.add_argument('--samples_yaml',
            type=str, required=True,
            help='path to the .yaml file containing sample info')
        sample_add_samples_yaml.set_defaults(
            func = self.add_samples_from_yaml_cmdline)

        # ADD SAMPLE
        sample_add = projects_subparser.add_parser('add_sample',
            description = help_desc['add_sample'],
            help = help_desc['add_sample'],
            parents = [self.__get_project_sample_parser(),
                     self.__get_sample_extra_arguments_parser(),
                     self.__get_barcode_flavor_species_parser(species_required=True),
                     self.__get_read_parser(reads_required=True)])
        sample_add.set_defaults(func=self.add_update_delete_sample_cmdline,
            action='add')

        # UPDATE SAMPLE
        sample_update = projects_subparser.add_parser('update_sample',
            description = help_desc['update_sample'],
            help = help_desc['update_sample'],
            parents = [self.__get_project_sample_parser(),
                     self.__get_sample_extra_arguments_parser(),
                     self.__get_barcode_flavor_species_parser(),
                     self.__get_read_parser()])
        sample_update.set_defaults(func=self.add_update_delete_sample_cmdline,
            action = 'update')

        # DELETE SAMPLE
        sample_remove = projects_subparser.add_parser('delete_sample',
            description = help_desc['delete_sample'],
            help = help_desc['delete_sample'],
            parents = [self.__get_project_sample_parser()])
        sample_remove.set_defaults(func=self.add_update_delete_sample_cmdline,
            action = 'delete')

        # MERGE SAMPLES
        sample_merge = projects_subparser.add_parser('merge_samples',
            description = help_desc['merge_samples'],
            help = help_desc['merge_samples'],
            parents=[self.__get_project_sample_lists_parser('perform the merge'),
                    self.__get_sample_extra_arguments_parser()])
        sample_merge.add_argument('--merged_project_id',
            required=True, type = str)
        sample_merge.add_argument('--merged_sample_id',
            required=True, type = str)
        sample_merge.set_defaults(func=self.merge_samples_cmdline)

        # LIST PROJECTS
        # always show these variables
        always_show = ['puck_id', 'species', 'investigator', 'sequencing_date',
            'experiment']
        remaining_options = [x for x in self.project_df_default_values.keys() \
            if x not in always_show]
        list_projects = projects_subparser.add_parser('list',
            description = help_desc['list'],
            help = help_desc['list'],
            parents=[self.__get_project_sample_lists_parser(
                'subset the data. If none provided, all will be displayed')])
        list_projects.add_argument('--variables',
            help = 'which extra variables to display per sample? ' +
                    f'{always_show} will always be shown.',
            choices = remaining_options,
            default=[],
            nargs='*')
        list_projects.set_defaults(func=self.list_projects_cmdline,
            always_show=always_show)
        
        return parser
