__version__ = "0.1"
__author__ = ["Nikos Karaiskos"]
__license__ = "GPL"

from spacemake.project_df import ProjectDF
from spacemake.config import ConfigFile
from spacemake.errors import SpacemakeError

include: 'scripts/snakemake_helper_functions.py'
include: 'variables.py'

project_df = ProjectDF(config['project_df'], ConfigFile.from_yaml('config.yaml'))

rule all:
    input:
        unpack(get_output_files)