import pytest

from spacemake.map_strategy import *
from spacemake.config import ConfigFile
from spacemake.project_df import ProjectDF
from spacemake.errors import *
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

def test_validation(test_root):
    config = ConfigFile.from_yaml((test_root / "config.yaml").as_posix())
    data = [
        ("flipped", "bowtie2:rRNA->STAR:genome", 'test_hsa', "rRNA:bowtie2->genome:STAR"),
        ("species_missing", "bowtie2:rRNA->STAR:genome", 'test_hs', "<class 'spacemake.errors.ConfigVariableNotFoundError'>"),
        ("with_cflavor", "rRNA:bowtie2@custom_index->genome:STAR@default", 'test_hsa', "rRNA:bowtie2@custom_index->genome:STAR@default"),
        ("unknown_cflavor", "rRNA:bowtie2@custom->genome:STAR@default", 'test_hsa', "<class 'spacemake.errors.ConfigVariableNotFoundError'>"),
        # ("flipped", "bowtie2:rRNA->STAR:genome", "rRNA:bowtie2->genome:STAR"),
    ]
    for name, mapstr, species, expect in data:
        try:
            res = validate_mapstr(mapstr, config=config, species=species)
        except (ValueError, ConfigVariableNotFoundError) as e:
            res = str(type(e))

        # print(f"test '{name}': {mapstr}-> {res} expect={expect} {expect == res}")
        assert res == expect


def test_mapstr(test_root):
    config = ConfigFile.from_yaml((test_root / "config.yaml").as_posix())
    data = [
        ("with_cflavor", "rRNA:bowtie2@custom_index->genome:STAR@default", None),
    ]
    for name, mapstr, expect in data:
        mr, lr = mapstr_to_targets(mapstr)
        assert mr[0].input_name == "uBAM"
        assert mr[-1].input_name == "rRNA.bowtie2"
        assert lr[0].link_src == "genome.STAR"
        assert lr[0].link_name == "final"
    

def test_get_mapped_BAM_output(test_root):
    config = ConfigFile.from_yaml((test_root / "config.yaml").as_posix())
    project_df = ProjectDF((test_root / "project_df.csv").as_posix(), config=config)

    out_files = get_mapped_BAM_output(project_df=project_df, config=config)
    print(out_files)
