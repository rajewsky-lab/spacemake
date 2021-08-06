import pandas as pd
import os
import errno
import yaml
import math
import argparse
import datetime
import re
from functools import reduce
from operator import getitem

class ConfigFile:
    line_separator = '-'*50+'\n'

    def __init__(self, file_path):
        self.file_path = file_path
        self.variables = yaml.load(open(file_path),
                    Loader=yaml.FullLoader)

    def dump(self):
        with open(self.file_path, 'w') as fo:
            fo.write(yaml.dump(self.variables))

    def set_file_path(self, file_path):
        self.file_path = file_path
        
    @property
    def puck_data(self):
        return self.variables['puck_data']

    def add_species_info(self, name, genome, annotation):
        if 'annotations' not in self.variables['knowledge']:
            self.variables['knowledge']['annotations'] = {}

        if 'genomes' not in self.variables['knowledge']:
            self.variables['knowledge']['genomes'] = {}

        if name in self.variables['knowledge']['annotations'].keys() and\
                name in self.variables['knowledge']['genomes'].keys():
           return False

        if os.path.isfile(annotation):
            self.variables['knowledge']['annotations'][name] = annotation
        else:
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), annotation)

        if os.path.isfile(genome):
            self.variables['knowledge']['genomes'][name] = genome
        else:
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), genome)

        return True

    def run_mode_exists(self, run_mode):
        return run_mode in self.variables['run_modes'].keys()

    def barcode_flavor_exists(self, barcode_flavor):
        return barcode_flavor in self.variables['knowledge']['barcode_flavor'].keys()
    
    def get_variable(self, path):
        return reduce(getitem, path, self.variables)  

    def list_variable(self, path):
        print(yaml.dump(self.get_variable(path)))

    def delete_run_mode_cmdline(self, args):
        name = args['name']
        if not name in self.variables['run_modes'].keys():
            print(f'run_mode: {name} does not exists, so it cannot be deleted')
            return 0
        del self.variables['run_modes'][name]
        self.dump()

    def list_run_modes_cmdline(self, args):
        self.list_variable(['run_modes'])
        
    def add_run_mode_cmdline(self, args):
        name = args['name']
        del args['name']
        
        # add and save new run mode
        if name in self.variables['run_modes']:
            msg = f'run_mode: {name} already exists, so it cannot be added.\n'
            msg += 'you can update this run_mode using the `spacemake config update_run_mode`'
            msg += ' command.\nor delete it using the `spacemake config delete_run_mode`' 
            msg += ' command.'

            print(msg)
            return 0

        msg = f'adding run_mode: {name}\n'
        msg += self.line_separator
        msg += 'variables:\n'
        msg += self.line_separator
        msg += yaml.dump(args, sort_keys=False)

        print(msg)

        self.variables['run_modes'][name] = args
        self.dump()

        print(f'run_mode: {name} added')

        return 1
    
    def update_run_mode_cmdline(self, args):
        name = args['name']
        del args['name']

        if not name in self.variables['run_modes']:
            msg = f'run_mode: {name} does not exists, so it cannot be updated.\n'
            msg += 'add this run_mode using the `spacemake config add_run_mode` command'
            print(msg)
            return 0

        msg = f'updating run_mode: {name}\n'
        msg += self.line_separator
        msg += 'variables:\n'
        msg += self.line_separator
        msg += yaml.dump(args, sort_keys=False)

        print(msg)

        self.variables['run_modes'][name].update(args)
        self.dump()

        print(f'run_mode: {name} updated')

        return 1

    def get_add_update_run_mode_parser(self):
        parser = argparse.ArgumentParser(
            description='add/update run_mode parent parser',
            add_help = False)
        parser.add_argument('name', type=str,
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
    
    def list_barcode_flavors_cmdline(self, args):
        self.list_variable(['knowledge', 'barcode_flavor'])

    def delete_barcode_flavor_cmdline(self, args):
        flavor_name = args['name']
        barcode_flavor = self.get_variable(['knowledge', 'barcode_flavor'])

        msg = self.line_separator
        msg += f'deleting {flavor_name} from barcode flavors\n'
        msg += self.line_separator

        if flavor_name not in barcode_flavor.keys():
            msg += 'barcode flavor with {flavor_name} do not exists'
            msg += ', so cant be deleted.\n'
            msg += 'aborting'
        else:
            flavor = barcode_flavor[flavor_name]
            msg += yaml.dump(flavor)
            msg += '\n' + self.line_separator
            del barcode_flavor[flavor_name]
            self.dump()
            msg += 'success!'

        print(msg)

    def add_barcode_flavor_cmdline(self, args):
        umi = args['umi']
        cell_barcode = args['cell_barcode']
        flavor_name = args['name']
        bam_tags = 'CR:{cell},XC:{cell},XM:{UMI},RG:{assigned}'

        # r(1|2) and then string slice
        to_match = r'r(1|2)(\[((?=-)-\d+|\d)*\:((?=-)-\d+|\d*)(\:((?=-)-\d+|\d*))*\])+$'

        msg = self.line_separator
        msg += f'adding {flavor_name} to barcode flavors\n'
        msg += self.line_separator

        barcode_flavor = self.get_variable(['knowledge', 'barcode_flavor'])

        if flavor_name in barcode_flavor.keys():
            msg += f'barcode flavor {flavor_name} already exists.\n'
            msg += 'to update you need to delete, and set it again.\n'
            msg += 'aborting'
            
            print(msg)
            return 1

        if re.match(to_match, umi) is None:
            msg += f'provided umi {umi} has the wrong structure.\n'
            msg += 'aborting'

            print(msg)
            return 1

        if re.match(to_match, cell_barcode) is None:
            msg += f'provided cell_barcode {cell_barcode} has the wrong structure.\n'
            msg += 'aborting'

            print(msg)
            return 1

        barcode_flavor[flavor_name] = {
            'UMI': umi,
            'bam_tags': bam_tags,
            'cell': cell_barcode}
        
        msg += yaml.dump(barcode_flavor[flavor_name])
        
        self.dump()

        msg += 'success!'
        print(msg)

        return 0

    def list_species_cmdline(self, args):
        # list annotations first
        msg = self.line_separator
        msg += 'annotations\n'
        msg += self.line_separator
        print(msg)
        self.list_variable(['knowledge', 'annotations'])

        # list genomes next
        msg = self.line_separator
        msg += 'genomes\n'
        msg += self.line_separator
        print(msg)
        self.list_variable(['knowledge', 'genomes'])

    def delete_species_cmdline(self, args):
        species_name = args['name']
        annotations = self.get_variable(['knowledge', 'annotations'])
        genomes = self.get_variable(['knowledge', 'genomes'])

        msg = self.line_separator
        msg += f'deleting {species_name} annotation and genome info.\n'
        msg += self.line_separator

        if species_name not in annotations.keys():
            msg += 'species {species_name} not in annotations.\n'
            msg += 'skipping.\n'
        else:
            del annotations[species_name]
            msg += 'deleted from annotations.\n'

        if species_name not in genomes.keys():
            msg += 'species {species_name} not in genomes.\n'
            msg += 'skipping.\n'
        else:
            del genomes[species_name]
            msg += 'deleted from genomes.\n'

        self.dump()

        msg += self.line_separator
        msg += 'success!'

        print(msg)

    def add_species_cmdline(self, args):
        species_name = args['name']
        annotation = args['annotation']
        genome = args['genome']

        msg = self.line_separator
        msg += f'adding genome, annotation of {species_name}\n'
        msg += self.line_separator

        try:
            added = self.add_species_info(species_name, genome, annotation)
            if added:
                msg += f'added {species_name}.\n'
                msg += f'genome: {genome}\n'
                msg += f'annotation: {annotation}\n'
                msg += self.line_separator
                msg += 'success!'
            else:
                msg += f'{species_name} already exists!\n'
                msg += 'aborting.'
        except FileNotFoundError as e:
            msg += f'Error: {e.filename} is not found!\n'
            msg += 'aborting.'

        self.dump()

        print(msg) 

    def get_subparsers(self, subparsers):
        parser_config = subparsers.add_parser('config', help = 'configure spacemake')
        parser_config_subparsers = parser_config.add_subparsers(help = 'config sub-command help')

        ## list run_modes ##
        parser_config_list_run_modes = parser_config_subparsers.add_parser('list_run_modes',
            help = 'list available run_modes')
        parser_config_list_run_modes.set_defaults(func=self.list_run_modes_cmdline)

        # add run_mode
        parser_config_add_run_mode = parser_config_subparsers.add_parser('add_run_mode',
            help = 'add a new run_mode',
            parents=[self.get_add_update_run_mode_parser()])
        parser_config_add_run_mode.set_defaults(func=self.add_run_mode_cmdline)

        # update run mode
        parser_config_update_run_mode = parser_config_subparsers.add_parser('update_run_mode',
            help = 'update run_mode',
            parents=[self.get_add_update_run_mode_parser()])
        parser_config_update_run_mode.set_defaults(func=self.update_run_mode_cmdline)

        parser_config_delete_run_mode = parser_config_subparsers.add_parser('delete_run_mode',
            help = 'delete a run_mode')
        parser_config_delete_run_mode.add_argument('name',
            type=str,
            help='run_mode to be deleted')
        parser_config_delete_run_mode.set_defaults(func=self.delete_run_mode_cmdline)

        # list barcode flavors
        parser_config_list_barcode_flavors = parser_config_subparsers\
            .add_parser('list_barcode_flavors',
                help = 'list barcode flavors and their settings')
        parser_config_list_barcode_flavors.set_defaults(
            func=self.list_barcode_flavors_cmdline)

        # delete barcode flavor
        parser_config_delete_barcode_flavor = parser_config_subparsers\
            .add_parser('delete_barcode_flavor',
                help = 'delete barcode flavor')
        parser_config_delete_barcode_flavor.add_argument('name',
            help = 'name of the barcode flavor to be deleted',
            type=str)
        parser_config_delete_barcode_flavor.set_defaults(
            func=self.delete_barcode_flavor_cmdline)

        # add barcode flavor
        parser_config_add_barcode_flavor = parser_config_subparsers\
            .add_parser('add_barcode_flavor',
                help = 'add a new barcode_flavor')
        parser_config_add_barcode_flavor.add_argument('name',
            help = 'name of the barcode flavor', type = str)
        parser_config_add_barcode_flavor.add_argument('umi',
            help = 'structure of UMI', type=str)
        parser_config_add_barcode_flavor.add_argument('cell_barcode',
            help = 'structure of CELL BARCODE', type=str)
        parser_config_add_barcode_flavor.set_defaults(
            func=self.add_barcode_flavor_cmdline)

        # list species settings
        parser_config_list_species = parser_config_subparsers\
            .add_parser('list_species',
                help = 'list all defined species and their genomes, annotations')
        parser_config_list_species.set_defaults(
            func=self.list_species_cmdline)

        # delete species
        parser_config_delete_species = parser_config_subparsers\
            .add_parser('delete_species',
                help = 'delete a species (genome and annotation) from configuration')
        parser_config_delete_species.add_argument('name',
            help = 'name of the species to be deleted',
            type=str)
        parser_config_delete_species.set_defaults(
            func=self.delete_species_cmdline)
        # add species
        parser_config_add_species = parser_config_subparsers\
            .add_parser('add_species',
                help = 'add a new species: genome (.fa) and annotation (.gtf) files')
        parser_config_add_species.add_argument('name',
            help = 'name of the species to be added',
            type = str)
        parser_config_add_species.add_argument('genome',
            help = 'path to the genome (.fa) file for the species to be added',
            type = str)
        parser_config_add_species.add_argument('annotation',
            help = 'path to the annotation (.gtf) file for the species to be added',
            type = str)
        parser_config_add_species.set_defaults(
            func=self.add_species_cmdline)

        return parser_config

class ProjectDF:
    # default values of the project dataframe columns
    line_separator = '-'*50+'\n'
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
        config: ConfigFile
    ):
        self.file_path = file_path
        self.config = config

        if os.path.isfile(file_path):
            self.df = pd.read_csv(file_path,
                index_col=['project_id', 'sample_id'],
                converters={'run_mode': eval})
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
    
    def add_sample_cmdline(self, args):
        project_id = args['project_id']
        sample_id = args['sample_id']

        msg = f'Adding ({project_id}, {sample_id})\n'
        msg += self.line_separator
        sample_added, sample = project_df.add_sample(**args)

        if sample_added :
            msg += 'SUCCESS: sample added successfully\n'
            msg += self.line_separator
            msg += sample.__str__()
            project_df.dump()
        else:
            msg += f'ERROR: sample with index ({project_id}, {sample_id})'
            msg +=' already exists. you need to update it first\n'
            msg +='use `spacemake projects update_sample` to update it'

        print(msg)

    def update_sample_cmdline(self, args):
        project_id = args['project_id']
        sample_id = args['sample_id']

        msg = f'Updating ({project_id}, {sample_id})\n'
        msg += self.line_separator
        sample_updated, sample = project_df.update_sample(**args)

        if sample_updated:
            msg += 'SUCCESS: sample updated successfully\n'
            msg += self.line_separator
            msg += sample.__str__()
            project_df.dump()
        else:
            msg += f'ERROR: sample with index ({project_id}, {sample_id})'
            msg +=' does not exist. you need to add it first\n'
            msg +='use `spacemake projects add_sample` to update it'

        print(msg)

    def get_add_sample_sheet_parser(self):
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

    def add_sample_sheet_cmdline(self, args):
        self.add_sample_sheet(args['sample_sheet'],
            args['basecalls_dir'])

        self.dump()

    def add_samples_from_yaml_cmdline(self, args):
        self.add_samples_from_yaml(args['samples_yaml'])

        self.dump()

    def set_species_cmdline(self, args):
        print(self.config.puck_data)
        projects = args['project_id']
        samples = args['sample_id']

        ix = self.df.query(
            'project_id in @projects or sample_id in @samples').index

        self.df.loc[ix, 'species'] = args['species_name']

        self.dump()

    def get_subparsers(self, subparsers):
        parser_project = subparsers.add_parser('projects', help ='manage projects and samples') 
        parser_sample_subparsers = parser_project.add_subparsers(help = 'sample sub-command help')

        # ADD SAMPLE SHEET
        parser_sample_add_sample_sheet = parser_sample_subparsers.add_parser('add_sample_sheet',
            help = 'add projects and samples from Illumina sample sheet',
            parents=[self.get_add_sample_sheet_parser()])
        parser_sample_add_sample_sheet.set_defaults(func=self.add_sample_sheet_cmdline)

        # ADD SAMPLES FROM YAML
        parser_sample_add_samples_yaml = parser_sample_subparsers.add_parser('add_samples_from_yaml',
            help = 'add several samples at once from a .yaml file')
        parser_sample_add_samples_yaml.add_argument('samples_yaml',
            type=str,
            help='path to the .yaml file containing sample info')
        parser_sample_add_samples_yaml.set_defaults(
            func = self.add_samples_from_yaml_cmdline)

        # ADD SAMPLE
        parser_sample_add = parser_sample_subparsers.add_parser('add_sample',
            help = 'add new sample',
            parents=[self.get_add_update_sample_parser(),
                     self.get_read_species_parser(reads_required=True)])
        parser_sample_add.set_defaults(func=self.add_sample_cmdline)

        # UPDATE SAMPLE
        parser_sample_add = parser_sample_subparsers.add_parser('update_sample',
            help = 'update_existing_sample',
            parents=[self.get_add_update_sample_parser(),
                     self.get_read_species_parser(reads_required=False)])
        parser_sample_add.set_defaults(func=self.update_sample_cmdline)

        # SET SPECIES
        parser_set_species = parser_sample_subparsers.add_parser('set_species',
            help = 'set species for one or multiple projects/samples')
        parser_set_species.add_argument('--species_name',
            help = 'name of the species to be added',
            type=str, required=True)
        parser_set_species.add_argument('--project_id',
            nargs='+',
            default=[],
            help = 'project id-s for which we should set this species')
        parser_set_species.add_argument('--sample_id',
            nargs='+',
            default=[],
            help = 'sample id-s for which we should set this species')
        parser_set_species.set_defaults(func=self.set_species_cmdline)
        
        return parser_project