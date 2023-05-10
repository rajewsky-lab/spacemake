import pytest
from spacemake.cmdline import *
import sys
import os

@pytest.fixture(scope="session")
def test_root(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("root")
    sm_path = os.path.dirname(__file__)
    # make a tmp-copy of the test_config.yaml
    def_config = os.path.join(sm_path, "../test_data/test_config.yaml")
    os.system(f"cp {def_config} {tmp / 'config.yaml'}")

    test_pdf =  os.path.join(sm_path, "../test_data/test_project_df.csv")
    os.system(f"cp {test_pdf} {tmp / 'project_df.csv'}")

    return tmp

def sm(*argc):
    sys.argv = ["spacemake",] + list(argc)
    cmdline()

def test_parsers(test_root):
    get_project_sample_parser()
    get_add_sample_sheet_parser()
    
    parser = get_sample_main_variables_parser()
    get_action_sample_parser(parser.add_subparsers(), "add", lambda *argc, **kw: None)

    get_sample_extra_info_parser()
    get_data_parser()
    get_run_parser()

    make_main_parser()


def test_init(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("disposable")
    os.chdir(tmp.as_posix())

    # just get the version
    sm("--version")

    # test the init parser
    sm("init")

def test_subcommands(test_root):
    # test everything after init
    os.chdir(test_root.as_posix())

    sm("projects", "list")
    sm("config", "list_adapter_flavors")

def test_species(test_root):
    os.chdir(test_root.as_posix())
    spacemake_dir = os.path.dirname(__file__) + '/../'

    sm(
        "config", "add_species",
        "--name=hsa_test",
        "--reference=genome",
        f"--annotation={spacemake_dir}/test_data/test_genome.fa.gz",
        f"--sequence={spacemake_dir}/test_data/test_genome.fa.gz",
    )

def test_sample(test_root):
    # test everything after init
    os.chdir(test_root.as_posix())
    spacemake_dir = os.path.dirname(__file__) + '/../'
    sm(
        "projects", "add_sample",
        "--project_id=test", "--sample_id=test1", 
        f'--R1={spacemake_dir}/test_data/reads_chr22_R1.fastq.gz', 
        f'--R2={spacemake_dir}/test_data/reads_chr22_R1.fastq.gz', 
        '--map_strategy="STAR:genome:final"',
        '--species=hsa_test'
    )

def test_run(test_root):
    # test everything after init
    # os.chdir(test_root.as_posix())
    # spacemake_dir = os.path.dirname(__file__) + '/../'
    sm("run", "-np", "--cores=8")

