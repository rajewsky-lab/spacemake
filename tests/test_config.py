import pytest
import spacemake.config as config
import os

@pytest.fixture(scope="session")
def test_yaml(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("config")
    sm_path = os.path.dirname(__file__)
    # make a tmp-copy of the test_config.yaml
    def_config = os.path.join(sm_path, "../test_data/test_config.yaml")
    test_config = tmp / "config.yaml"
    os.system(f"cp {def_config} {test_config}")

    return test_config

def test_parsers():
    config.get_puck_parser()
    config.get_run_mode_parser()
    config.get_species_parser()
    config.get_barcode_flavor_parser()
    config.get_adapter_parser()
    config.get_adapter_flavor_parser()


def test_default(test_yaml):
    cfg = config.ConfigFile.from_yaml(test_yaml.as_posix())
    cfg.add_variable("adapters", name="test_adapter", seq="ACGTACGTACGTAAAGGG")
    cfg.add_variable("barcode_flavors", name="test_bc_flavor", umi="r1[12:20] + r2[:8]", cell="r1[:12]")
    cfg.add_variable("adapter_flavors", name="test_adapter_flavor", cut_left=["SMART:max_errors=0.2:min_overlap=3"], cut_right=["Q:min_base_qual=20", "polyG:max_errors=0.2:min_overlap=3", "polyA:max_errors=0.2:min_overlap=3"])
    
    sm_path = os.path.dirname(__file__)
    test_data = os.path.join(sm_path, "../test_data/")
    cfg.add_variable("species", name="test_species", reference="genome", annotation=test_data + "/gencode.v38.chr22.gtf", sequence=test_data + "/hg38_chr22.fa", STAR_flags="")
    cfg.add_variable("species", name="test_species", reference="rRNA", sequence=test_data + "/rRNA_hsa.fa.gz", BT2_flags="--very-sensitive-local")







