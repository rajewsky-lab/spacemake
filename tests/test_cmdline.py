import pytest
from spacemake.cmdline import *
import sys

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


def test_parsers(test_root):
    get_project_sample_parser()
    get_add_sample_sheet_parser()
    
    parser = get_sample_main_variables_parser()
    get_action_sample_parser(parser.add_subparsers(), "add", lambda *argc, **kw: None)

    get_sample_extra_info_parser()
    get_data_parser()
    get_run_parser()


def test_init(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("disposable")
    os.chdir(tmp.as_posix())

    # just get the version
    sys.argv = ["spacemake", "--version"]
    # args = make_main_parser().parse_args()
    # args.func()
    cmdline()

    # test the init parser
    sys.argv = ["spacemake", "init"]
    cmdline()

def test_subcommands(test_root):
    # test everything after init
    os.chdir(test_root.as_posix())
    sys.argv = ["spacemake", "--version"]
    cmdline()

    sys.argv = ["spacemake", "projects", "list"]
    cmdline()

    sys.argv = ["spacemake", "config", "list_adapter_flavors"]
    cmdline()

    
