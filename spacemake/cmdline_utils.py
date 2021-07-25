import pandas as pd
import os
import yaml
import math
import argparse
import datetime

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

    def add_species_info(self, species_info):
        for species, data_type in species_info.items():
            if not species in self.variables['knowledge']:
                self.variables['knowledge'][species] = {}

            for data_type_name, data_type_value in data_type.items():
                if data_type_name in self.variables['knowledge'][species].keys():
                    print(f'{species} already has {data_type_name}, omitting.')
                self.variables['knowledge'][species][data_type_name] = data_type_value

    def run_mode_exists(self, run_mode):
        return run_mode in self.variables['run_modes'].keys()

    def barcode_flavor_exists(self, barcode_flavor):
        return barcode_flavor in self.variables['knowledge']['barcode_flavor'].keys()

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
        cf = cls(args['config_file_path'])

        run_modes = cf.variables['run_modes']

        print(yaml.dump(run_modes))

        return run_modes

    @classmethod
    def add_run_mode(cls, args):
        cf = cls(args['config_file_path'])

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
        cf = cls(args['config_file_path'])

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

    @staticmethod 
    def get_add_update_run_mode_parser():
        parser = argparse.ArgumentParser(
            description='add/update run_mode parent parser',
            add_help = False)
        parser.add_argument('--name', type=str, required=True,
            help='name of the run_mode to be added')
        parser.add_argument('--parent_run_mode', type=str,
            help='Name of the parent run_mode. All run_modes will fall back to \'default\'')
        parser.add_argument('--umi_cutoff', type=int, nargs='+',
            help='umi_cutoff for this run_mode.' +\
                 'the automated analysis will be run with these cutoffs')
        parser.add_argument(
            '--clean_dge',
            required=False,
            action=argparse.BooleanOptionalAction,
            help='if set, the DGE will be cleaned of barcodes which overlap with primers')
        parser.add_argument(
            '--detect_tissue',
            required=False,
            action=argparse.BooleanOptionalAction,
            help='By default only beads having at least umi_cutoff UMI counts are analysed '+\
                 'during the automated analysis, all other beads are filtered out. If this ' +\
                 'parameter is set, contiguous islands within umi_cutoff passing beads will '+\
                 'also be included in the analysis')
        parser.add_argument(
            '--plot_bead_size',
            help='The bead size to be used when plotting the beads in 2D, during the '+\
                 'automated report generation. Defaults to 1.')
        parser.add_argument(
            '--polyA_adapter_trimming',
            required=False,
            action=argparse.BooleanOptionalAction,
            help='If set, reads will have polyA stretches and adapter sequence overlaps trimmed '+\
                 'BEFORE mapping.')
        parser.add_argument(
            '--count_intronic_reads',
            required=False,
            action=argparse.BooleanOptionalAction,
            help='If set, INTRONIC reads will also be countsed (apart from UTR and CDS)')
        parser.add_argument(
            '--count_mm_reads',
            required=False,
            action=argparse.BooleanOptionalAction,
            help='If True, multi-mappers will also be counted. For every multimapper only reads which '+\
                 'have one unique read mapped to a CDS or UTR region will be counted')

        return parser

    @staticmethod
    def add_config_subparsers(subparsers, config_path):
        parser_config = subparsers.add_parser('config', help = 'configure spacemake')
        parser_config_subparsers = parser_config.add_subparsers(help = 'config sub-command help')

        ## list run_modes ##
        parser_config_list_run_modes = parser_config_subparsers.add_parser('list_run_modes',
                help = 'list available run_modes')
        parser_config_list_run_modes.set_defaults(func=ConfigFile.list_run_modes,
                config_file_path = config_path)

        # add run_mode
        parser_config_add_run_mode = parser_config_subparsers.add_parser('add_run_mode',
                help = 'add a new run_mode', parents=[ConfigFile.get_add_update_run_mode_parser()])
        parser_config_add_run_mode.set_defaults(func=ConfigFile.add_run_mode,
                config_file_path = config_path)

        # update run mode
        parser_config_update_run_mode = parser_config_subparsers.add_parser('update_run_mode',
                help = 'update run_mode',
                parents=[ConfigFile.get_add_update_run_mode_parser()])
        parser_config_update_run_mode.set_defaults(func=ConfigFile.update_run_mode,
                config_file_path = config_path)

        parser_config_delete_run_mode = parser_config_subparsers.add_parser('delete_run_mode',
            help = 'delete a run_mode')
        parser_config_delete_run_mode.add_argument('--name',
            required=True,
            type=str,
            help='run_mode to be deleted')
        parser_config_delete_run_mode.set_defaults(func=ConfigFile.delete_run_mode)
        
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
        "is_merged":False}
    
    def __init__(
        self,
        file_path,
        puck_data = {
            'barcode_file': 'barcode_file.csv',
            'root': 'puck_data'
        }
    ):
        self.file_path = file_path

        if os.path.isfile(file_path):
            self.df = pd.read_csv(file_path,
                index_col=['project_id', 'sample_id'])
        else:
            index = pd.MultiIndex(
                names=['project_id', 'sample_id'],
                levels=[[],[]],
                codes=[[],[]]) 
            self.df = pd.DataFrame(columns = self.project_df_default_values.keys(),
                index=index)

        self.puck_data = puck_data

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

        puck_dir = find_dir(puck_id, self.puck_data['root'])
        path = None

        if puck_dir is not None:
            # puck dir exists, look for barcode file pattern
            path = os.path.join(puck_dir, self.puck_data["barcode_file"])

            return get_barcode_file(path)
        else:
            return self.project_df_default_values['puck_barcode_file']

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

        if self.sample_exists(project_id, sample_id):
            print(f'sample with {ix} already exists in ProjectDF')
            print(f'updating')
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
            return (False, None)  
        else:
            return (True,
                    self.__add_update_sample(project_id, sample_id,
                    return_series=True, **kwargs))

    def update_sample(self, project_id, sample_id, **kwargs):
        if self.sample_exists(project_id, sample_id):
            return (True, self.__add_update_sample(project_id, sample_id,
                    return_series=True, **kwargs))
        else:
            return (False, None)

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

        #project_df = df_assign_merge_samples(project_df)
    
    @staticmethod
    def get_add_update_sample_parser():
        parser = argparse.ArgumentParser(
            add_help=False)

        parser.add_argument('--project_id',
            type=str, required=True,
            help = 'project_id of the sample to be added/updated')

        parser.add_argument('--sample_id',
            type=str, required=True,
            help = 'sample_id of the sample to be added/updated')

        parser.add_argument('--puck_id',
            type=str,
            help = 'puck_id of the sample to be added/update')

        parser.add_argument('--puck_barcode_file',
            type=str,
            help = 'the path to the file contining (x,y) positions of the barcodes')

        parser.add_argument('--investigator',
            type=str,
            help = 'add the name of the investigator(s) responsible for this sample')

        parser.add_argument('--sequencing_date',
            type=datetime.date.fromisoformat,
            help = 'sequencing date of the sample. format: YYYY-MM-DD')

        parser.add_argument('--run_mode',
            type=str,
            nargs='+',
            help = 'run_mode names for this sample. the sample will be processed using the provided run_modes')

        parser.add_argument('--barcode_flavor',
            type=str,
            help = 'barcode flavor for this sample')

        return parser

    @staticmethod
    def get_read_species_parser(reads_required):
        parser = argparse.ArgumentParser(
            add_help = False)

        parser.add_argument('--R1',
            type=str,
            help = '.fastq.gz file path to R1 reads',
            required=reads_required)

        parser.add_argument('--R2',
            type=str,
            help = '.fastq.gz file path to R2 reads',
            required=reads_required)

        parser.add_argument('--species',
            type=str,
            help = 'add the name of the species for this sample',
            required=reads_required)


        return parser
    
    @classmethod
    def add_sample_cmdline(cls, args):
        project_id = args['project_id']
        sample_id = args['sample_id']

        msg = f'Adding ({project_id}, {sample_id})\n'
        msg += '----------------------------------\n'
        project_df = cls(args['project_df_file'])
        sample_added, sample = project_df.add_sample(**args)

        if sample_added :
            msg += 'SUCCESS: sample added successfully\n'
            msg += '----------------------------------\n'
            msg += sample.__str__()
            project_df.dump()
        else:
            msg += f'ERROR: sample with index ({project_id}, {sample_id})'
            msg +=' already exists. you need to update it first\n'
            msg +='use `spacemake projects update_sample` to update it'

        print(msg)

    @classmethod
    def update_sample_cmdline(cls, args):
        project_id = args['project_id']
        sample_id = args['sample_id']

        msg = f'Updating ({project_id}, {sample_id})\n'
        msg += '----------------------------------\n'
        project_df = cls(args['project_df_file'])
        sample_updated, sample = project_df.update_sample(**args)

        if sample_updated:
            msg += 'SUCCESS: sample updated successfully\n'
            msg += '----------------------------------\n'
            msg += sample.__str__()
            project_df.dump()
        else:
            msg += f'ERROR: sample with index ({project_id}, {sample_id})'
            msg +=' does not exist. you need to add it first\n'
            msg +='use `spacemake projects add_sample` to update it'

        print(msg)

    @staticmethod
    def get_add_sample_sheet_parser():
        parser = argparse.ArgumentParser(
            description = 'add a new sample sheet to the samples',
            add_help=False)

        parser.add_argument('--sample_sheet', type = str,
            help = 'the path to the Illumina sample sheet',
            required=True)
        parser.add_argument('--basecalls_dir', type = str,
            help = 'path to the basecalls directory',
            required=True)

        return parser

    @classmethod
    def add_sample_sheet_cmdline(cls, args):
        project_df = cls(args['project_df_file'])
        project_df.add_sample_sheet(args['sample_sheet'],
            args['basecalls_dir'])

        project_df.dump()

    @classmethod
    def add_samples_from_yaml_cmdline(cls, args):
        project_df = cls(args['project_df_file'])

        project_df.add_samples_from_yaml(args['samples_yaml'])

        project_df.dump()

