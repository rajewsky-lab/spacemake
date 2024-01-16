import pytest
import os

spacemake_dir = os.path.abspath(os.path.dirname(__file__) + '/../')

def sm(*argc, expect_fail=False):
    # construct the desired cmdline
    import sys
    sys.argv = ["spacemake",] + list(argc)
    
    # ensure that no ConfigFile and ProjectDF instances 
    # are retained from previous tests
    import spacemake.config
    import spacemake.project_df
    spacemake.config.__global_config = None
    spacemake.project_df.__global_ProjectDF = None

    # execute spacemake cmdline code
    from spacemake.cmdline import cmdline
    res = cmdline()
    if expect_fail:
        assert isinstance(res, Exception) == True
    else:
        assert isinstance(res, Exception) == False

    return res

@pytest.fixture
def tmp_root(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("root_blank")

    return tmp


@pytest.fixture
def initialized_root(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("root_initialized")
    os.chdir(tmp.as_posix())

    # just get the version
    sm("--version")

    # test the init parser
    sm("init", "--dropseq_tools", "/data/rajewsky/shared_bins/Drop-seq_tools-2.5.1/")

    return tmp


@pytest.fixture
def with_species(initialized_root):
    os.chdir(initialized_root.as_posix())
    # # test old way
    # sm(
    #     "config", "add_species",
    #     "--name=hsa_test",
    #     f"--genome={spacemake_dir}/test_data/test_genome.fa.gz",
    #     f"--annotation={spacemake_dir}/test_data/test_annotation.gtf.gz",
    # )
    # test new way
    sm(
        "config", "add_species",
        "--name=test_hsa",
        "--reference=genome",
        f"--sequence={spacemake_dir}/test_data/test_genome.fa.gz",
        f"--annotation={spacemake_dir}/test_data/test_annotation.gtf.gz",
    )
    # add a second reference
    sm(
        "config", "add_species",
        "--name=test_hsa",
        "--reference=rRNA",
        f"--sequence={spacemake_dir}/test_data/rRNA_hsa.fa.gz",
    )
    # add a third reference
    sm(
        "config", "add_species",
        "--name=test_hsa",
        "--reference=miRNA",
        f"--sequence={spacemake_dir}/test_data/mirgenedb.hsa.mature.fa.gz",
    )
    return initialized_root


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
