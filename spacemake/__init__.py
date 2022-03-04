__version__ = 1.0

from . import preprocess as pp
from . import spatial as sp
from .config import ConfigFile
from .project_df import ProjectDF
import argparse

def bead_image_generation_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)

    parser.add_argument('--project_id', type=str,
        required=True)

    parser.add_argument('--sample_id', type=str,
        required=True)

    parser.add_argument('--run_mode', type=str,
        required=True)

    return parser

class SpaceMake:
    def __init__(self, root):
        self.root = root
        self.config = ConfigFile.from_yaml(f'{root}/config.yaml')
        self.project_df = ProjectDF(f'{root}/project_df.csv', config=self.config)
    
    def load_processed_adata(self,
        project,
        sample,
        run_mode_name,
        umi_cutoff 
    ):
        run_mode = self.config.get_run_mode(run_mode_name)
        
        adata = sc.read(f'{self.root}/projects/{project}/processed_data/{sample}/'+
                f'illumina/complete_data/automated_analysis/{run_mode_name}/' +
                f'umi_cutoff_{umi_cutoff}/results.h5ad')
        
        adata.uns['run_mode_variables'] = run_mode.variables
        adata.uns['puck_variables'] = self.project_df.get_puck_variables(
            project_id = project,
            sample_id = sample)
        
        return adata
    
    def load_raw_spatial_adata(
        self,
        project,
        sample,
        run_mode_name
    ):
        run_mode = self.config.get_run_mode(run_mode_name)

        dge_type = ""
        dge_cleaned = ""
        polyA_adapter_trimmed = ""
        mm_included = ""

        if run_mode.variables["polyA_adapter_trimming"]:
            polyA_adapter_trimmed = ".polyA_adapter_trimmed"

        if run_mode.variables["count_intronic_reads"]:
            dge_type = ".all"
        else:
            dge_type = ".exon"

        if run_mode.variables["count_mm_reads"]:
            mm_included = ".mm_included"

        if run_mode.variables["clean_dge"]:
            dge_cleaned = ".cleaned"

        adata = sc.read(f'{self.root}/projects/{project}/processed_data/{sample}/'+
                f'illumina/complete_data/dge/dge{dge_type}{dge_cleaned}'+
                f'{polyA_adapter_trimmed}{mm_included}.spatial_beads.h5ad')
        
        adata.uns['run_mode_variables'] = run_mode.variables
        adata.uns['puck_variables'] = self.project_df.get_puck_variables(
            project_id = project,
            sample_id = sample)
        
        return adata
