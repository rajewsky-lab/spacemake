import pytest
import os


@pytest.fixture
def tmp_root(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("root_blank")

    return tmp

@pytest.fixture
def configured_root(tmp_path_factory):
    tmp_root = tmp_path_factory.mktemp("root_preconfigured")

    sm_path = os.path.dirname(__file__)
    # make a tmp-copy of the test_config.yaml
    def_config = os.path.join(sm_path, "../test_data/test_config.yaml")
    os.system(f"cp {def_config} {tmp_root / 'config.yaml'}")

    test_pdf =  os.path.join(sm_path, "../test_data/test_project_df.csv")
    os.system(f"cp {test_pdf} {tmp_root / 'project_df.csv'}")

    return tmp_root
