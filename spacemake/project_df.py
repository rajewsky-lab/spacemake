import pandas as pd
import os
import yaml
import math
import argparse
import datetime
import re
import logging

from spacemake.errors import *
from spacemake.config import ConfigFile
from spacemake.util import message_aggregation, assert_file

logger_name = 'spacemake.project_df'

def get_project_sample_parser(allow_multiple=False, prepend='', help_extra=''):
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False)

    project_argument = 'project_id'
    sample_argument = 'sample_id'
    required = True
    nargs = None
    default = None

    if allow_multiple:
        project_argument = f'{project_argument}_list'
        sample_argument = f'{sample_argument}_list'
        required = False
        nargs = '*'
        default = []

    parser.add_argument(f'--{prepend}{project_argument}',
        type=str, required=required, nargs=nargs, default=default,
        help = f'{project_argument}. {help_extra}')

    parser.add_argument(f'--{prepend}{sample_argument}',
        type=str, required=required, nargs=nargs, default=default,
        help = f'{sample_argument}. {help_extra}')

    return parser

def get_add_sample_sheet_parser():
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

def get_sample_extra_arguments_parser(species_required=False,
        reads_required=False):
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False)

    parser.add_argument('--R1',
        type=str,
        help = '.fastq.gz file path to R1 reads',
        required=reads_required)

    parser.add_argument('--R2',
        type=str,
        help = '.fastq.gz file path to R2 reads',
        required=reads_required)

    parser.add_argument('--barcode_flavor',
        type=str,
        help = 'barcode flavor for this sample')

    parser.add_argument('--species',
        type=str,
        help = 'add the name of the species for this sample',
        required=species_required)

    parser.add_argument('--puck', type=str,
        help = 'name of the puck for this sample. if puck contains a \n'+
            '`barcodes` path to a coordinate file, those coordinates\n' +
            ' will be used when processing this sample. if \n' +
            ' not provided, a default puck will be used with \n'+
            'width_um=3000, spot_diameter_um=10')

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

    return parser

def get_set_remove_variable_subparsers(parent_parser, 
    variable_name, func, prepend = '', allow_multiple=False):
    if allow_multiple:
        nargs = '+'
    else:
        nargs=None

    def get_action_parser(action):
        action_parser = parent_parser.add_parser(f'{action}_{variable_name}',
            parents=[get_project_sample_parser(allow_multiple=True)],
            description = f'{action} {variable_name} for several projects/samples',
            help = f'{action} {variable_name} for several projects/samples')
        action_parser.add_argument(f'--{variable_name}',
            type = str, required=True,
            nargs = nargs)
        action_parser.set_defaults(
            func=func, variable=variable_name,
            action = action)

        return action_parser

    set_parser = get_action_parser('set')

    if allow_multiple:
        set_parser.add_argument('--keep_old', action='store_true')

        remove_parser = get_action_parser('remove')

def get_action_sample_parser(parent_parser, action, func):
    if action not in ['add', 'update', 'delete', 'merge']:
        raise ValueError(f'Invalid action: {action}')

    if action == 'merge':
        parser_name = 'merge_samples'
        msg = 'merge samples'
        parents = [get_project_sample_parser(prepend='merged_'),
                   get_project_sample_parser(allow_multiple=True)]
    else:
        parser_name = f'{action}_sample'
        msg = f'{action} a sample'
        parents = [get_project_sample_parser()]


    if action == 'add':
        parents.append(get_sample_extra_arguments_parser(
            species_required = True, reads_required = True))
    elif action == 'update' or action == 'merge':
        parents.append(get_sample_extra_arguments_parser())

    sample_parser = parent_parser.add_parser(parser_name,
        description = msg, 
        help = msg,
        parents = parents)
    sample_parser.set_defaults(
        func=func,
        action=action)

def setup_project_parser(pdf, attach_to):
    parser = attach_to.add_parser('projects',
            help='manage projects and samples',
            description ='Using one of the subcommands specified below, it is possible to' +\
                    ' add/update/remove projects and their settings') 
    subparsers = parser.add_subparsers()

    help_desc = {
        'add_sample_sheet' : 'add projects and samples from Illumina sample sheet',
        'add_samples_from_yaml' : 'add several samples at once from a .yaml file',
        'merge_samples': 'merge several samples into one. ' +\
                'samples need to have the same species',
        'list':'list several project(s) and sample(s) and their settings'
    }
    
    # ADD/UPDATE/DELETE/MERGE sample(s)
    # get parsers for adding, deleting and updating a sample
    for action in ['add', 'update', 'delete']:
        get_action_sample_parser(subparsers, action,
            # attach the function to the parser
            func = lambda args: add_update_delete_sample_cmdline(pdf, args))

    # get parser for merging
    get_action_sample_parser(subparsers, 'merge',
        func = lambda args: merge_samples_cmdline(pdf, args))

    # LIST PROJECTS
    # always show these variables
    always_show = ['puck_id', 'species', 'investigator', 'sequencing_date',
        'experiment']
    remaining_options = [x for x in pdf.project_df_default_values.keys() \
        if x not in always_show]
    list_projects = subparsers.add_parser('list',
        description = help_desc['list'],
        help = help_desc['list'],
        parents=[get_project_sample_parser(allow_multiple=True,
            help_extra='subset the data. If none provided, all will be displayed')])
    list_projects.add_argument('--variables',
        help = 'which extra variables to display per sample? ' +
                f'{always_show} will always be shown.',
        choices = remaining_options,
        default=[],
        nargs='*')
    list_projects.set_defaults(func=lambda args: list_projects_cmdline(pdf, args),
        always_show=always_show)

    # ADD SAMPLE SHEET
    sample_add_sample_sheet = subparsers.add_parser('add_sample_sheet',
        description = help_desc['add_sample_sheet'],
        help = help_desc['add_sample_sheet'],
        parents=[get_add_sample_sheet_parser()])
    sample_add_sample_sheet.set_defaults(
        func=lambda args: add_sample_sheet_cmdline(pdf, args))

    # ADD SAMPLES FROM YAML
    sample_add_samples_yaml = subparsers.add_parser('add_samples_from_yaml',
        description = help_desc['add_samples_from_yaml'],
        help = help_desc['add_samples_from_yaml'])
    sample_add_samples_yaml.add_argument('--samples_yaml',
        type=str, required=True,
        help='path to the .yaml file containing sample info')
    sample_add_samples_yaml.set_defaults(
        func = lambda args: add_samples_from_yaml_cmdline(pdf, args))

    # get set/remove parser for each main variable
    # this will add parser for:
    # setting/removing run_mode
    # setting species
    # setting puck
    # setting barcode_flavor
    for main_var_sg, main_var_type in ConfigFile.main_variable_sg2type.items():
        allow_multiple = False
        if isinstance(main_var_type, str) and main_var_type.endswith('_list'):
            allow_multiple = True

        get_set_remove_variable_subparsers(subparsers, 
            variable_name = main_var_sg,
            func=lambda args: set_remove_variable_cmdline(pdf, args),
            allow_multiple=allow_multiple)

    return parser

@message_aggregation(logger_name)
def add_update_delete_sample_cmdline(pdf, args):
    action = args['action']

    # remove the action from args
    del args['action']

    if action == 'add' or action =='update':
        func = lambda **kwargs: pdf.add_update_sample(action=action, **kwargs)
    elif action == 'delete':
        func = pdf.delete_sample

    sample = func(**args)
    pdf.dump()

@message_aggregation(logger_name)
def add_sample_sheet_cmdline(pdf, args):
    pdf.add_sample_sheet(args['sample_sheet'],
        args['basecalls_dir'])

    pdf.dump()

@message_aggregation(logger_name)
def add_samples_from_yaml_cmdline(pdf, args):
    pdf.add_samples_from_yaml(args['samples_yaml'])

    pdf.dump()

@message_aggregation(logger_name)
def set_remove_variable_cmdline(pdf, args):
    variable_name = args['variable']
    action = args['action']
    variable_key = args[variable_name]

    pdf.set_remove_variable(
        variable_name = variable_name,
        variable_key = variable_key,
        action = action,
        project_id_list = args['project_id_list'],
        sample_id_list = args['sample_id_list'],
        keep_old = args.get('keep_old', False))

    pdf.dump()

@message_aggregation(logger_name)
def merge_samples_cmdline(pdf, args):
    pdf.merge_samples(**args)

    pdf.dump()

@message_aggregation(logger_name)
def list_projects_cmdline(pdf, args):
    projects = args['project_id_list']
    samples = args['sample_id_list']
    variables = args['always_show'] + args['variables']

    df = pdf.df
    logger = logging.getLogger(logger_name)

    if projects != [] or samples != []:
        df = df.query('project_id in @projects or sample_id in @samples')
        logger.inof(f'listing projects: {projects} and samples: {samples}')
    else:
        logger.info('listin all projects and samples')

    logger.info(f'variables used: {variables}')
    
    # print the table
    logger.info(df.loc[:, variables].__str__())

class ProjectDF:
    # default values of the project dataframe columns
    project_df_default_values = {
        "puck_id": "no_optical_puck",
        "sample_sheet": None,
        "species": None,
        "demux_barcode_mismatch": 1,
        "demux_dir": None,
        "basecalls_dir": None,
        "R1": None,
        "R2": None,
        "investigator": "unknown",
        "sequencing_date": "unknown",
        "experiment": "unknown",
        "puck_barcode_file": None,
        "run_mode": ["default"],
        "barcode_flavor": "default",
        "is_merged":False,
        "merged_from":[],
        "puck":None}
    
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
                converters={'run_mode': eval, 'merged_from': eval},
                na_values = ['None', 'none'])
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

        self.logger = logging.getLogger(logger_name)

    def compute_max_barcode_mismatch(self, indices):
        """computes the maximum number of mismatches allowed for demultiplexing based
        on the indices present in the sample sheet."""
        num_samples = len(indices)

        if num_samples == 1:
            return 4
        else:
            max_mismatch = 3
            for i in range(num_samples - 1):
                for j in range(i + 1, num_samples):
                    hd = self.hamming_distance(indices[i], indices[j])
                    max_mismatch = min(max_mismatch, math.ceil(hd / 2) - 1)

        return max_mismatch

    def hamming_distance(self, string1, string2):
        return sum(c1 != c2 for c1, c2 in zip(string1, string2))
    
    def find_barcode_file(self, puck_id):
        # first find directory of puck file

        # return none or the path of the file
        def get_barcode_file(path):
            if os.path.isfile(path):
                return path

            return None

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

    def is_spatial(self, project_id ,sample_id):
        puck_barcode_file = self.get_metadata('puck_barcode_file',
            project_id = project_id,
            sample_id = sample_id)

        puck_name = self.get_metadata('puck',
            project_id = project_id,
            sample_id = sample_id)

        puck = self.config.get_puck(puck_name, return_empty=True)
        
        if puck_barcode_file is not None or puck.has_barcodes:
            return True
        else:
            return False

    def get_puck_variables(self, project_id, sample_id, return_empty=False):
        puck_name = self.get_metadata('puck',
            project_id = project_id,
            sample_id = sample_id)

        return self.config.get_puck(puck_name, return_empty=return_empty).variables

    def get_metadata(self, field, project_id=None, sample_id=None, **kwargs):
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
            investigator = None
            sequencing_date = None

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
        df["demux_barcode_mismatch"] = self.compute_max_barcode_mismatch(df["index"])
        df["sample_sheet"] = sample_sheet_path
        df["demux_dir"] = df["sample_sheet"].str.split("/").str[-1].str.split(".").str[0]
        df["puck_barcode_file"] = df.puck_id.apply(self.find_barcode_file)
        df.set_index(['project_id', 'sample_id'], inplace=True)

        for ix, row in df.iterrows():
            self.add_update_sample(project_id=ix[0], sample_id=ix[1], **row.to_dict())

    def sample_exists(self, project_id = None, sample_id = None):
        if project_id is None or sample_id is None:
            raise Exception(f'you need to provide a sample_id and project_id to check if sample exists')
        else:
            ix = (project_id, sample_id) 

            return ix in self.df.index

    def add_update_sample(self, action=None, project_id = None, sample_id = None,
            R1=None, R2=None, sample_sheet=None, basecalls_dir=None, is_merged=False,
            return_series = False,
            **kwargs):
        """
        adds or updates a sample with a given project_id and sample_id
        """
        print(kwargs)
        ix = (project_id, sample_id)
        sample_exists = self.sample_exists(*ix)

        if action is None and sample_exists:
            action = 'update'
        elif action is None and not sample_exists:
            action = 'add'
        elif action == 'add' and sample_exists:
            raise SampleAlreadyExistsError(ix)
        elif action == 'update' and not sample_exists:
            raise ProjectSampleNotFoundError('(project_id, sample_id)',ix)

        if action == 'add':
            self.logger.info(f'Adding sample {ix}')
        elif action == 'update':
            self.logger.info(f'Updating sample {ix}')
        else:
            raise ValueError(f'Unknown action {action}')
        
        # check variables
        if action == 'add' and (R1 is None or R2 is None):
            self.logger.info('R1 or R2 not provided, trying basecalls_dir and'+
                ' sample_sheet')
            if basecalls_dir is None or sample_sheet is None:
                raise ValueError('Neither R1 and R2 nor basecalls_dir and ' +
                    'sample_sheet weren\'t provided. you need to provide the'
                    ' first two or the second two')

        # assert files first
        assert_file(R1, default_value = None, extension = '.fastq.gz')

        assert_file(R2, default_value = None, extension = '.fastq.gz')

        is_spatial = assert_file(
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

        # if everything correct, add or update
        # first populate kwargs
        kwargs['R1'] = R1
        kwargs['R2'] = R2
        kwargs['sample_sheet'] = sample_sheet
        kwargs['basecalls_dir'] = basecalls_dir
        kwargs['is_merged'] = is_merged

        if sample_exists:
            new_project = self.df.loc[ix].copy()
            # pd.Series.update will only update values which are not None
            new_project.update(pd.Series(kwargs))
            self.df.loc[ix] = new_project
        else:
            # if sample is spatial, and puck not provided, assign 'default'
            if is_spatial and 'puck' not in kwargs:
                kwargs['puck'] = 'default'

            new_project = pd.Series(self.project_df_default_values)
            new_project.name = ix
            new_project.update(kwargs)

            self.df = self.df.append(new_project)

        if return_series:
            return (ix, new_project)

    def delete_sample(self, project_id, sample_id):
        if self.sample_exists(project_id, sample_id):
            self.logger.info(f'Deleting sample: {ix}, with the following' +
                    f' with the following variables:\n{element}')

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

        if demux_projects is not None:
            # if we have projects in the config file
            # get the samples
            for ip in demux_projects:
                self.add_sample_sheet(ip['sample_sheet'], ip['basecalls_dir'])

        # add additional samples from config.yaml, which have already been demultiplexed.
        for project in config['additional_projects']:
            self.add_update_sample(**project)

    def get_ix_from_project_sample_list(self,
        project_id_list = [], sample_id_list = []):
        # raise error if both lists are empty
        if project_id_list == [] and sample_id_list == []:
            raise NoProjectSampleProvidedError()

        # of only one provided use that, if both use intersection
        if project_id_list == []:
            ix = self.df.query(
                'sample_id in @sample_id_list').index
        elif sample_id_list == []:
            ix = self.df.query(
                'project_id in @project_id_list').index
        else:
            ix = self.df.query(
                'project_id in @project_id_list and sample_id in @sample_id_list').index

        return ix

    def set_variable(self, ix, variable_name, variable_key, keep_old=False):
        variable_name_pl = self.config.main_variables_sg2pl[variable_name]
        self.config.assert_main_variable(variable_name_pl)
        self.config.assert_variable(variable_name_pl, variable_key)
        
        # get current value
        i_variable = self.df.at[ix, variable_name]

        # if list, we either append or not
        if isinstance(i_variable, list):
            if keep_old:
                # keep the old list as well
                i_variable = list(set(i_variable + variable_key))
            else:
                i_variable = list(set(variable_key))

        # if we do not keep the old, simply set
        else:
            i_variable = variable_key

        if i_variable == [] or i_variable is None:
            raise EmptyConfigVariableError(variable_name)

        self.df.at[ix, variable_name] = i_variable


    def remove_variable(self, ix, variable_name, variable_key):
        variable_name_pl = self.config.main_variables_sg2pl[variable_name]
        self.config.assert_main_variable(variable_name_pl)
        self.config.assert_variable(variable_name_pl, variable_key)

        if not isinstance(variable_key, list):
            raise ValueError('variable_key has to be a list')

        i_variable = self.df.at[ix, variable_name]

        i_variable = [val for val in i_variable if val not in variable_key]
        if i_variable == [] or i_variable is None:
            raise EmptyConfigVariableError(variable_name)

        self.df.at[ix, variable_name] = i_variable

    def set_remove_variable(self, variable_name, variable_key,
            action, project_id_list = [], sample_id_list = [], keep_old = False):
        self.assert_projects_samples_exist(project_id_list, sample_id_list)

        variable_name_pl = self.config.main_variables_sg2pl[variable_name]
        self.config.assert_main_variable(variable_name_pl)

        ix = self.get_ix_from_project_sample_list(
            project_id_list = project_id_list,
            sample_id_list = sample_id_list)

        if isinstance(variable_key, list):
            for var in variable_key:
                self.config.assert_variable(variable_name_pl, var)
        else:
            self.config.assert_variable(variable_name_pl, variable_key)
        
        for i, row in self.df.loc[ix, :].iterrows():
            # add/remove variable. if it already exists dont do anything
            if action == 'set':
                self.set_variable(i, variable_name, variable_key, keep_old)
            elif action == 'remove':
                self.remove_variable(i, variable_name, variable_key)

        return ix.to_list()


    def assert_projects_samples_exist(self,
            project_id_list = [], sample_id_list = []):
        for project in project_id_list:
            if project not in self.df.index.get_level_values('project_id'):
                raise ProjectSampleNotFoundError('project_id', project)

        for sample in sample_id_list:
            if sample not in self.df.index.get_level_values('sample_id'):
                raise ProjectSampleNotFoundError('sample_id', sample)
        
    def merge_samples(self,
        merged_project_id,
        merged_sample_id,
        project_id_list = [],
        sample_id_list = [],
        **kwargs
    ):
        # check if projects and samples with these IDs exist
        self.assert_projects_samples_exist(project_id_list, sample_id_list)

        ix = self.get_ix_from_project_sample_list(
            project_id_list = project_id_list,
            sample_id_list = sample_id_list)

        consistent_variables = list(self.config.main_variables_sg2pl.keys())
        consistent_variables.remove('run_mode')

        ix_list = ix.to_list()

        if ix_list == []:
            raise ProjectSampleNotFoundError(
                '(project_id_list, sample_id_list)',
                (project_id_list, sample_id_list))

        # check for variable inconsistency
        # raise error if variable different between samples
        for variable in consistent_variables:
            if variable in kwargs:
                # if variable provided from command line, skip
                continue

            variable_val = self.df.loc[ix, variable].to_list()
            
            if len(set(variable_val)) > 1:
                raise InconsistentVariablesDuringMerge(
                    variable_name = variable,
                    variable_value = variable_val,
                    ix = ix)
            else:
                # attach the deduced, consisten variable
                kwargs[variable] = variable_val[0]

        # get puck_barcode_file
        if 'puck_barcode_file' not in kwargs:
            pbf_default = self.project_df_default_values['puck_barcode_file']
            pbf_list = self.df.loc[ix, 'puck_barcode_file'].to_list()
            # filter out default values
            pbf_list = [x for x in pbf_list if x != pbf_default]

            # remove duplicates
            pbf_list = list(set(pbf_list))            

            if pbf_list == []:
                # if all values are default
                kwargs['puck_barcode_file'] = pbf_default
            elif len(pbf_list) == 1:
                kwargs['puck_barcode_file'] = pbf_list[0]
            else:
                raise InconsistentVariablesDuringMerge(
                    variable_name = 'puck_barcode_file',
                    variable_value = pbf_list,
                    ix = ix.to_list())


        variables_to_deduce = [
            'investigator',
            'experiment',
            'sequencing_date'
        ]

        for variable in variables_to_deduce:
            if variable not in kwargs:
                kwargs[variable] = ';'.join(self.df.loc[ix, variable].unique())

            
        # if no run_mode provided, overwrite with user defined one
        if 'run_mode' not in kwargs.keys():
            run_mode_lists = [set(run_mode)
                for run_mode in self.df.loc[ix].run_mode.to_list()
            ]

            # join run modes from parent samples
            if len(run_mode_lists) == 1:
                run_mode = run_mode_lists[0]
            else:
                run_mode = run_mode_lists[0].intersection(*run_mode_lists[1:])

            # create a list from the set intersection
            run_mode = list(run_mode)

            # if there are no common elements, throw an error
            if len(run_mode) == 0:
                raise InconsistentVariablesDuringMerge(
                    variable_name = 'run_mode',
                    variable_value = run_mode_lists,
                    ix = ix_list)
            
            # finally add run mode to arguments
            kwargs['run_mode'] = run_mode
        
        sample_added = self.add_update_sample(
            action = 'add',
            project_id = merged_project_id,
            sample_id = merged_sample_id,
            is_merged=True,
            merged_from = ix.to_list(),
            **kwargs)

        return (sample_added, ix)