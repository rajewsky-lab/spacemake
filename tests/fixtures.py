import pytest
import os

spacemake_dir = os.path.abspath(os.path.dirname(__file__) + "/../")
print("SPACEMAKE_DIR", spacemake_dir)


def sm(*argc, expect_fail=False):
    # construct the desired cmdline
    import sys

    sys.argv = [
        "spacemake",
    ] + list(argc)

    # ensure that no ConfigFile and ProjectDF instances
    # are retained from previous tests
    import spacemake.config
    import spacemake.project_df

    spacemake.config.__global_config = None
    spacemake.project_df.__global_ProjectDF = None

    # execute spacemake cmdline code
    from spacemake.cmdline import cmdline

    res = cmdline()
    # print("res", res)
    if expect_fail:
        assert isinstance(res, Exception) == True
    else:
        assert isinstance(res, Exception) == False

    return res


def _init():
    # just get the version
    sm("--version")

    # test the init parser
    sm("init", "--dropseq-tools", "/data/rajewsky/shared_bins/Drop-seq_tools-2.5.1/")


def _add_species():
    sm(
        "config",
        "add-species",
        "--name=test_hsa",
        "--reference=genome",
        f"--sequence={spacemake_dir}/test_data/test_genome.fa.gz",
        f"--annotation={spacemake_dir}/test_data/test_genome.gtf.gz",
    )
    # add a second reference
    sm(
        "config",
        "add-species",
        "--name=test_hsa",
        "--reference=rRNA",
        f"--sequence={spacemake_dir}/test_data/rRNA_hsa.fa.gz",
    )
    # add a third reference
    sm(
        "config",
        "add-species",
        "--name=test_hsa",
        "--reference=miRNA",
        f"--sequence={spacemake_dir}/test_data/mirgenedb.hsa.mature.fa.gz",
    )
    # pretend we have mouse as well
    # TODO: place some actual mouse genome and/or phiX genomes in test-data repository
    sm(
        "config",
        "add-species",
        "--name=mouse",
        "--reference=genome",
        f"--sequence={spacemake_dir}/test_data/test_genome.fa.gz",
        f"--annotation={spacemake_dir}/test_data/test_genome.gtf.gz",
    )
    sm(
        "config",
        "add-species",
        "--name=mouse",
        "--reference=phiX",
        f"--sequence={spacemake_dir}/test_data/test_genome.fa.gz",
        f"--annotation={spacemake_dir}/test_data/test_genome.gtf.gz",
    )
    sm(
        "config",
        "add-species",
        "--name=mouse",
        "--reference=rRNA",
        f"--sequence={spacemake_dir}/test_data/rRNA_hsa.fa.gz",
    )


@pytest.fixture
def tmp_root(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("root_blank")

    return tmp


@pytest.fixture
def initialized_root(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("root_initialized")
    os.chdir(tmp.as_posix())

    _init()
    return tmp


@pytest.fixture
def with_species(initialized_root):
    os.chdir(initialized_root.as_posix())
    # # test old way
    # sm(
    #     "config", "add_species",
    #     "--name=hsa_test",
    #     f"--genome={spacemake_dir}/test_data/test_genome.fa.gz",
    #     f"--annotation={spacemake_dir}/test_data/test_genome.gtf.gz",
    # )
    # test new way
    _add_species()
    return initialized_root


@pytest.fixture
def configured_root(tmp_path_factory):
    tmp_root = tmp_path_factory.mktemp("root_preconfigured")

    # make a tmp-copy of the test_config.yaml
    def_config = os.path.join(spacemake_dir, "test_data/test_config.yaml")
    os.system(f"cp {def_config} {tmp_root / 'config.yaml'}")

    test_pdf = os.path.join(spacemake_dir, "test_data/test_project_df.csv")
    open(f"{tmp_root / 'project_df.csv'}", "w").write(
        open(test_pdf, "r").read().format(spacemake_dir=spacemake_dir)
    )
    # os.system(f"cp {test_pdf} {tmp_root / 'project_df.csv'}")

    return tmp_root


@pytest.fixture(scope="session")
def with_tile_test_data(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("root_tile_test")
    os.chdir(tmp.as_posix())
    _init()
    _add_species()
    print(
        "return code",
        os.system(
            "wget https://bimsbstatic.mdc-berlin.de/rajewsky/spacemake-test-data/spacemake_tile_test_data.tar.gz -O /dev/stdout | tar -xz"
        ),
    )
    print(os.listdir("."))

    return tmp
