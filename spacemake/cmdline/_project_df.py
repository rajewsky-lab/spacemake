import pandas as pd
import os
import yaml
import math
import argparse
import datetime
import re

from ..errors import *
from ._config import ConfigFile
from ._util import assert_file, LINE_SEPARATOR

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
        "puck":"none"}
    
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

    def is_spatial(self, project_id ,sample_id):
        puck_barcode_file = self.get_metadata('puck_barcode_file',
            project_id = project_id,
            sample_id = sample_id)

        puck_name = self.get_metadata('puck',
            project_id = project_id,
            sample_id = sample_id)

        puck = self.config.get_puck(puck_name, return_empty=True)
        
        if puck_barcode_file != 'none' or puck.has_barcodes:
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

        if self.sample_exists(project_id, sample_id):
            new_project = self.df.loc[ix].copy()
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

        if demux_projects is not None:
            # if we have projects in the config file
            # get the samples
            for ip in demux_projects:
                self.add_sample_sheet(ip['sample_sheet'], ip['basecalls_dir'])

        # add additional samples from config.yaml, which have already been demultiplexed.
        for project in config['additional_projects']:
            self.__add_update_sample(**project)

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
            help = 'name of the puck for this sample. if puck contains a \n'+
                '`barcodes` path to a coordinate file, those coordinates\n' +
                ' will be used when processing this sample. if \n' +
                ' not provided, a default puck will be used with \n'+
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

    def __set_variable(self, ix, variable_name, variable_key, keep_old=False):
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


    def __remove_variable(self, ix, variable_name, variable_key):
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

    def __set_remove_variable(self, variable_name, variable_key,
            action, project_id_list = [], sample_id_list = [], keep_old = False):
        self.__assert_projects_samples_exist(project_id_list, sample_id_list)

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
                self.__set_variable(i, variable_name, variable_key, keep_old)
            elif action == 'remove':
                self.__remove_variable(i, variable_name, variable_key)

        return ix.to_list()

    def set_remove_variable_cmdline(self, args):
        variable_name = args['variable']
        action = args['action']
        variable_key = args[variable_name]


        if action == 'set':
            succ_msg = 'set'
            action_msg = 'Setting'
        elif action == 'remove':
            succ_msg = 'removed'
            action_msg = 'Removing'

        msg =''

        try:
            set_indices = self.__set_remove_variable(
                variable_name = variable_name,
                variable_key = variable_key,
                action = action,
                project_id_list = args['project_id_list'],
                sample_id_list = args['sample_id_list'],
                keep_old = args.get('keep_old', False))

            msg += f'{action_msg} {variable_name}={variable_key} for '
            msg += f'the following (project_id, sample_id)'
            msg += f' samples:\n{set_indices}\n'
            msg += LINE_SEPARATOR
            msg += f'SUCCESS: {variable_name}={variable_key} {succ_msg} succesfully.\n'

            self.dump()
        except SpacemakeError as e:
            msg += str(e)
        finally:
            print(msg)

    def __set_remove_variable_subparsers(self, parent_parser, 
        variable_name, allow_multiple=False):
        if allow_multiple:
            nargs = '+'
        else:
            nargs=None

        def get_action_parser(action):
            action_parser = parent_parser.add_parser(f'{action}_{variable_name}',
                parents=[self.__get_project_sample_lists_parser()],
                description = f'{action} {variable_name}',
                help = f'{action} {variable_name}')
            action_parser.add_argument(f'--{variable_name}',
                type = str, required=True,
                nargs = nargs)
            action_parser.set_defaults(
                func=self.set_remove_variable_cmdline, variable=variable_name,
                action = action)

            return action_parser

        set_parser = get_action_parser('set')

        if allow_multiple:
            set_parser.add_argument('--keep_old', action='store_true')

            remove_parser = get_action_parser('remove')

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
        
    def merge_samples(self,
        merged_project_id,
        merged_sample_id,
        project_id_list = [],
        sample_id_list = [],
        **kwargs
    ):
        # check if projects and samples with these IDs exist
        self.__assert_projects_samples_exist(project_id_list, sample_id_list)

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

        main_config_vars = list(self.config.main_variables_sg2pl.keys())

        msg = f'\n\nthe following variables if not provided, will be checked'
        msg += f' for consistency between the merged sampels: \n'
        msg += f'{main_config_vars}'


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
                    self.__get_sample_extra_arguments_parser()],
            formatter_class=argparse.RawTextHelpFormatter)
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

        self.__set_remove_variable_subparsers(projects_subparser, 
            'species', allow_multiple=False)

        self.__set_remove_variable_subparsers(projects_subparser, 
            'run_mode', allow_multiple=True)
        
        self.__set_remove_variable_subparsers(projects_subparser, 
            'puck', allow_multiple=False)
        
        self.__set_remove_variable_subparsers(projects_subparser, 
            'barcode_flavor', allow_multiple=False)
        
        return parser
