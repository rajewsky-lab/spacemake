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
from spacemake.errors import FileWrongExtensionError, RunModeNotFoundError, \
    BarcodeFlavorNotFoundError, NoProjectSampleProvidedError, SpeciesNotFoundError, \
    ProjectSampleNotFoundError, SampleAlreadyExistsError, UnmatchingSpeciesDuringMerge

LINE_SEPARATOR = '-'*50+'\n'

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
        
    @property
    def puck_data(self):
        return self.variables['puck_data']

    def __get_knowledge_var(self, var_name):
        if var_name not in self.variables['knowledge']:
            self.variables['knowledge'][var_name] = {}

        return self.variables['knowledge'][var_name]

    @property
    def available_species(self):
        pass

    def __add_genome_annotation(self, var_name, species_name, file_path):
        var_values = self.__get_knowledge_var(var_name)

        if var_name == 'genomes' or var_name == 'rRNA_genomes':
            file_extension = '.fa'
        elif var_name == 'annotations':
            file_extension = '.gtf'

        assert_file(file_path, extension = file_extension)

        self.variables['knowledge'][var_name][species_name] = file_path

    def species_exists(self, species_name):
        if species_name in self.__get_knowledge_var('annotations').keys() \
                and species_name in self.__get_knowledge_var('genomes').keys():
            return True
        else:
            return False

    def add_species_info(self, name, genome, annotation,
            rRNA_genome=None):
        if self.species_exists(name):
           return False

        self.__add_genome_annotation('annotations', name, annotation)
        self.__add_genome_annotation('genomes', name, genome)

        if rRNA_genome is not None :
            self.__add_genome_annotation('rRNA_genomes', name, rRNA_genome)

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
        # set the name and delete from dictionary
        name = args['name']
        del args['name']
        
        msg = f'Adding run_mode: {name}\n'
        msg += LINE_SEPARATOR
        
        # add and save new run mode
        if name in self.variables['run_modes']:
            msg += f'ERROR: run_mode: {name} already exists, so it cannot be added.\n'
            msg += 'you can update this run_mode using the `spacemake config update_run_mode`'
            msg += ' command.\nor delete it using the `spacemake config delete_run_mode`' 
            msg += ' command.'
        else:
            self.variables['run_modes'][name] = args
            self.dump()
            msg += 'SUCCESS: run_mode added sucessfully.\n'
            msg += f'variables added for run_mode: {name}:\n'
            msg += LINE_SEPARATOR
            msg += yaml.dump(args, sort_keys=False)

        print(msg)
    
    def update_run_mode_cmdline(self, args):
        name = args['name']
        del args['name']

        if not name in self.variables['run_modes']:
            msg = f'run_mode: {name} does not exists, so it cannot be updated.\n'
            msg += 'add this run_mode using the `spacemake config add_run_mode` command'
            print(msg)
            return 0

        msg = f'updating run_mode: {name}\n'
        msg += LINE_SEPARATOR
        msg += 'variables:\n'
        msg += LINE_SEPARATOR
        msg += yaml.dump(args, sort_keys=False)

        print(msg)

        self.variables['run_modes'][name].update(args)
        self.dump()

        print(f'run_mode: {name} updated')

        return 1

    def __get_add_update_run_mode_parser(self):
        parser = argparse.ArgumentParser(
            allow_abbrev=False,
            description='add/update run_mode parent parser',
            add_help = False)
        parser.add_argument('--name', type=str,
            help='name of the run_mode to be added')
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
                 'automated report generation. Defaults to 1.', type=float)
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

        msg = LINE_SEPARATOR
        msg += f'deleting {flavor_name} from barcode flavors\n'
        msg += LINE_SEPARATOR

        if flavor_name not in barcode_flavor.keys():
            msg += 'barcode flavor with {flavor_name} do not exists'
            msg += ', so cant be deleted.\n'
            msg += 'aborting'
        else:
            flavor = barcode_flavor[flavor_name]
            msg += yaml.dump(flavor)
            msg += '\n' + LINE_SEPARATOR
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

        msg = LINE_SEPARATOR
        msg += f'adding {flavor_name} to barcode flavors\n'
        msg += LINE_SEPARATOR

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
        msg = LINE_SEPARATOR
        msg += 'annotations\n'
        msg += LINE_SEPARATOR
        print(msg)
        self.list_variable(['knowledge', 'annotations'])

        # list genomes next
        msg = LINE_SEPARATOR
        msg += 'genomes\n'
        msg += LINE_SEPARATOR
        print(msg)
        self.list_variable(['knowledge', 'genomes'])

    def delete_species_cmdline(self, args):
        species_name = args['name']
        annotations = self.get_variable(['knowledge', 'annotations'])
        genomes = self.get_variable(['knowledge', 'genomes'])

        msg = LINE_SEPARATOR
        msg += f'deleting {species_name} annotation and genome info.\n'
        msg += LINE_SEPARATOR

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

        msg += LINE_SEPARATOR
        msg += 'success!'

        print(msg)

    def add_species_cmdline(self, args):
        species_name = args['name']
        annotation = args['annotation']
        genome = args['genome']
        rRNA_genome = args.get('rRNA_genome', None)

        msg = LINE_SEPARATOR
        msg += f'adding genome, annotation of {species_name}\n'
        msg += LINE_SEPARATOR

        try:
            added = self.add_species_info(species_name, genome, annotation,
                rRNA_genome)
            if added:
                msg += f'added {species_name}.\n'
                msg += f'genome: {genome}\n'
                msg += f'annotation: {annotation}\n'
                if rRNA_genome is not None:
                    msg += f'rRNA_genome: {rRNA_genome}\n'
                msg += LINE_SEPARATOR
                msg += 'success!'
            else:
                msg += f'{species_name} already exists!\n'
                msg += 'aborting.'

            self.dump()
        except FileNotFoundError as e:
            msg += LINE_SEPARATOR
            msg += str(e)
        finally:
            print(msg) 

    def get_subparsers(self, subparsers):
        parser_config = subparsers.add_parser('config', help = 'configure spacemake')
        parser_config_subparsers = parser_config.add_subparsers(help = 'config sub-command help')

        ## list run_modes ##
        parser_config_list_run_modes = parser_config_subparsers.add_parser('list_run_modes',
            description = 'list available run_modes')
        parser_config_list_run_modes.set_defaults(func=self.list_run_modes_cmdline)

        # add run_mode
        parser_config_add_run_mode = parser_config_subparsers.add_parser('add_run_mode',
            description = 'add a new run_mode',
            parents=[self.__get_add_update_run_mode_parser()])
        parser_config_add_run_mode.set_defaults(func=self.add_run_mode_cmdline)

        # update run mode
        parser_config_update_run_mode = parser_config_subparsers.add_parser('update_run_mode',
            description = 'update run_mode',
            parents=[self.__get_add_update_run_mode_parser()])
        parser_config_update_run_mode.set_defaults(func=self.update_run_mode_cmdline)

        parser_config_delete_run_mode = parser_config_subparsers.add_parser('delete_run_mode',
            description = 'delete a run_mode')
        parser_config_delete_run_mode.add_argument('--name',
            type=str,
            help='run_mode to be deleted', required=True)
        parser_config_delete_run_mode.set_defaults(func=self.delete_run_mode_cmdline)

        # list barcode flavors
        parser_config_list_barcode_flavors = parser_config_subparsers\
            .add_parser('list_barcode_flavors',
                description = 'list barcode flavors and their settings')
        parser_config_list_barcode_flavors.set_defaults(
            func=self.list_barcode_flavors_cmdline)

        # delete barcode flavor
        parser_config_delete_barcode_flavor = parser_config_subparsers\
            .add_parser('delete_barcode_flavor',
                description = 'delete barcode flavor')
        parser_config_delete_barcode_flavor.add_argument('--name',
            help = 'name of the barcode flavor to be deleted',
            type=str, required=True)
        parser_config_delete_barcode_flavor.set_defaults(
            func=self.delete_barcode_flavor_cmdline)

        # add barcode flavor
        parser_config_add_barcode_flavor = parser_config_subparsers\
            .add_parser('add_barcode_flavor',
                description = 'add a new barcode_flavor')
        parser_config_add_barcode_flavor.add_argument('--name',
            help = 'name of the barcode flavor', type = str,
            required = True)
        parser_config_add_barcode_flavor.add_argument('--umi',
            help = 'structure of UMI', type=str, required=True)
        parser_config_add_barcode_flavor.add_argument('--cell_barcode',
            help = 'structure of CELL BARCODE', type=str, required=True)
        parser_config_add_barcode_flavor.set_defaults(
            func=self.add_barcode_flavor_cmdline)

        # list species settings
        parser_config_list_species = parser_config_subparsers\
            .add_parser('list_species',
                description = 'list all defined species and their genomes, annotations')
        parser_config_list_species.set_defaults(
            func=self.list_species_cmdline)

        # delete species
        parser_config_delete_species = parser_config_subparsers\
            .add_parser('delete_species',
                description = 'delete a species (genome and annotation) from configuration')
        parser_config_delete_species.add_argument('--name',
            help = 'name of the species to be deleted',
            type=str, required=True)
        parser_config_delete_species.set_defaults(
            func=self.delete_species_cmdline)
        # add species
        parser_config_add_species = parser_config_subparsers\
            .add_parser('add_species',
                description = 'add a new species: genome (.fa) and annotation (.gtf) files')
        parser_config_add_species.add_argument('--name',
            help = 'name of the species to be added',
            type = str, required=True)
        parser_config_add_species.add_argument('--genome',
            help = 'path to the genome (.fa) file for the species to be added',
            type = str, required=True)
        parser_config_add_species.add_argument('--annotation',
            help = 'path to the annotation (.gtf) file for the species to be added',
            type = str, required=True)
        parser_config_add_species.add_argument('--rRNA_genome',
            help = 'path to the ribosomal-RNA genome (.fa) file for the species to be added',
            default=None,
            type = str)
        parser_config_add_species.set_defaults(
            func=self.add_species_cmdline)

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
        "merged_from":[]}
    
    def __init__(
        self,
        file_path,
        config: ConfigFile = None
    ):
        self.file_path = file_path
        self.config = config

        if os.path.isfile(file_path):
            self.df = pd.read_csv(file_path,
                index_col=['project_id', 'sample_id'],
                converters={'run_mode': eval, 'merged_from': eval})
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
            if not self.config.run_mode_exists(run_mode):
                raise RunModeNotFoundError(run_mode)

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

    def remove_sample(self, project_id, sample_id):
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

        msg = ''

        if action == 'add':
            msg += 'Adding'
            succ_msg = 'added'
            func = self.add_sample
        elif action == 'update':
            msg += 'Updatig'
            succ_msg = 'updated'
            func = self.update_sample
        elif action == 'delete':
            msg += 'Deleting'
            succ_msg = 'deleted'
            func = self.remove_sample

        msg += f' ({project_id}, {sample_id})\n'

        try:
            sample = func(**args)

            msg += f'SUCCESS: sample {succ_msg} successfully.\n'
            msg += LINE_SEPARATOR
            msg += str(sample)
            self.dump()
        except (RunModeNotFoundError, FileNotFoundError,
                SampleAlreadyExistsError, ProjectSampleNotFoundError,
                FileWrongExtensionError, BarcodeFlavorNotFoundError) as e:
            msg += LINE_SEPARATOR
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

        if self.config.species_exists(species):
            self.df.loc[ix, 'species'] = species

            return ix.to_list()
        else:
            raise SpeciesNotFoundError(species)

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
        except (NoProjectSampleProvidedError, SpeciesNotFoundError,
                ProjectSampleNotFoundError) as e:
            msg += str(e)
        finally:
            print(msg)


    def __add_run_mode(self, ix, run_mode):
        i_run_mode = self.df.at[ix, 'run_mode']
        i_run_mode.append(run_mode)
        self.df.at[ix, 'run_mode'] = list(set(i_run_mode))

    def __remove_run_mode(self, ix, run_mode):
        i_run_mode = [rm for rm in self.df.at[ix, 'run_mode'] if rm != run_mode]
        self.df.at[ix, 'run_mode'] = i_run_mode

    def __set_remove_run_mode(self, run_modes, action, projects = [], samples = []):
        # raise error if both lists are empty
        if projects == [] and samples == []:
            raise NoProjectSampleProvidedError()

        self.__assert_projects_samples_exist(projects, samples)

        ix = self.df.query(
            'project_id in @projects or sample_id in @samples').index

        for run_mode in run_modes:
            if self.config.run_mode_exists(run_mode):
                for i, row in self.df.loc[ix, :].iterrows():
                    # add/remove run mode. if it already exists dont do anything
                    if action == 'set':
                        self.__add_run_mode(i, run_mode)
                    elif action == 'remove':
                        self.__remove_run_mode(i, run_mode)
            else:
                raise RunModeNotFoundError(run_mode)

        return ix.to_list()

    def set_remove_run_mode_cmdline(self, args):
        print(args)
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
        except (NoProjectSampleProvidedError, ProjectSampleNotFoundError, 
                RunModeNotFoundError) as e:

            msg += str(e)
        finally:
            print(msg)

    def list_projects_cmdline(self, args):
        projects = args['project_id_list']
        samples = args['sample_id_list']
        variables = args['variables']

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

        if self.config.barcode_flavor_exists(barcode_flavor):
            ix = self.df.query(
                'project_id in @project_id_list or sample_id in @sample_id_list').index

            if ix is None:
                raise NoProjectSampleProvidedError()

            self.df.loc[ix, 'barcode_flavor'] = barcode_flavor

            return ix.to_list()
        else:
            raise BarcodeFlavorNotFoundError(barcode_flavor)

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

        except (BarcodeFlavorNotFoundError, NoProjectSampleProvidedError,
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

        species = self.df.loc[ix].species.to_list()

        if ix.to_list() == []:
            raise ProjectSampleNotFoundError(
                '(project_id_list, sample_id_list)',
                (project_id_list, sample_id_list))

        if len(set(species)) > 1:
            raise UnmatchingSpeciesDuringMerge()
        else:
            species = species[0]

        
        fields_to_deduce = [
            'investigator',
            'experiment',
            'sequencing_date'
        ]

        for field in fields_to_deduce:
            if field not in kwargs.keys():
                kwargs[field] = ';'.join(self.df.loc[ix, field].unique())

            
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
                raise Exception('No run modes shared between merged samples')
            
            # finally add run mode to arguments
            kwargs['run_mode'] = run_mode
        
        sample_added = self.add_sample(
            project_id = merged_project_id,
            sample_id = merged_sample_id,
            species = species,
            is_merged=True,
            merged_from = ix.to_list(),
            **kwargs)

        return (sample_added, ix)

    def merge_samples_cmdline(self, args):
        msg = ''
        print(args)
        try:
            sample, set_indices = self.merge_samples(**args)

            set_indices = set_indices.to_list()

            msg += f'Merging samples {set_indices}.\n'
            msg += LINE_SEPARATOR
            msg += 'SUCCESS: samples merged successfully.\n'
            msg += LINE_SEPARATOR
            msg += str(sample)

            self.dump()
        except (NoProjectSampleProvidedError, SpeciesNotFoundError,
                ProjectSampleNotFoundError) as e:
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
            'set_species':'set species for one or multiple projects/samples',
            'set_barcode_flavor':'set barcode_flavor for one or multiple projects/samples',
            'set_run_mode':'set run mode(s) for one or multiple projects/samples',
            'remove_run_mode':'remove run mode(s) for one or multiple projects/samples',
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

        # SET SPECIES
        set_species = projects_subparser.add_parser('set_species',
            description = help_desc['set_species'],
            help = help_desc['set_species'],
            parents=[self.__get_project_sample_lists_parser('set the species')])

        set_species.add_argument('--species_name',
            help = 'name of the species to be set',
            type=str, required=True)
        set_species.set_defaults(func=self.set_species_cmdline)

        # SET BARCODE FLAVOR
        set_barcode_flavor = projects_subparser.add_parser('set_barcode_flavor',
            description = help_desc['set_barcode_flavor'],
            help = help_desc['set_barcode_flavor'],
            parents=[self.__get_project_sample_lists_parser('set barcode_flavor')])
        set_barcode_flavor.add_argument('--barcode_flavor',
            help = 'name of the barcode_flavor to be set',
            type=str, required=True)
        set_barcode_flavor.set_defaults(func=self.set_barcode_flavor_cmdline)

        # SET RUN MODE
        set_run_mode = projects_subparser.add_parser('set_run_mode',
            description = help_desc['set_run_mode'],
            help = help_desc['set_run_mode'],
            parents=[self.__get_project_sample_lists_parser('set run mode(s)')])
        set_run_mode.add_argument('--run_mode',
            help = 'name of the run mode(s) to be set',
            nargs = '+',
            type=str, required=True)
        set_run_mode.set_defaults(func=self.set_remove_run_mode_cmdline,
            action = 'set')
        
        # REMOVE RUN MODE
        remove_run_mode = projects_subparser.add_parser('remove_run_mode',
            description = help_desc['remove_run_mode'],
            help = help_desc['remove_run_mode'],
            parents=[self.__get_project_sample_lists_parser('remove run mode(s)')])
        remove_run_mode.add_argument('--run_mode',
            help = 'name of the run mode(s) to be removed',
            nargs = '+',
            type=str, required=True)
        remove_run_mode.set_defaults(func=self.set_remove_run_mode_cmdline,
            action = 'remove')

        # LIST PROJECTS
        list_projects = projects_subparser.add_parser('list',
            description = help_desc['list'],
            help = help_desc['list'],
            parents=[self.__get_project_sample_lists_parser(
                'subset the data. If none provided, all will be displayed')])
        list_projects.add_argument('--variables',
            help = 'which variables to display per sample? If no variables are provided,' +\
                    'puck_id, species, investigator, sequencing_date and experiment are shown.',
            choices = self.project_df_default_values.keys(),
            default = ['puck_id', 'species', 'investigator', 'sequencing_date', 'experiment'],
            nargs='*')
        list_projects.set_defaults(func=self.list_projects_cmdline)
        
        return parser
