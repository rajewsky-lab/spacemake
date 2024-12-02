import pytest

from spacemake.map_strategy import *
from spacemake.config import ConfigFile
from spacemake.project_df import ProjectDF
from spacemake.errors import *
import os

from fixtures import configured_root, tmp_root, sm, spacemake_dir


def test_validation(configured_root):
    config = ConfigFile.from_yaml((configured_root / "config.yaml").as_posix())
    data = [
        (
            "flipped",
            "rRNA:bowtie2->genome:STAR",
            "test_hsa",
            "bowtie2:rRNA->STAR:genome",
        ),
        (
            "species_missing",
            "bowtie2:rRNA->STAR:genome",
            "test_hs",
            "<class 'spacemake.errors.ConfigVariableNotFoundError'>",
        ),
        (
            "with_cflavor",
            "bowtie2@custom_index:rRNA->STAR@default:genome",
            "test_hsa",
            "bowtie2@custom_index:rRNA->STAR@default:genome",
        ),
        (
            "unknown_cflavor",
            "rRNA:bowtie2@custom->genome:STAR@default",
            "test_hsa",
            "<class 'spacemake.errors.ConfigVariableNotFoundError'>",
        ),
        # ("flipped", "bowtie2:rRNA->STAR:genome", "rRNA:bowtie2->genome:STAR"),
    ]
    for name, mapstr, species, expect in data:
        # print(f"running test {name}")
        try:
            res = validate_mapstr(mapstr, config=config, species=species)
        except (ValueError, ConfigVariableNotFoundError) as e:
            res = str(type(e))

        print(f"test '{name}': {mapstr}-> {res} expect={expect} {expect == res}")
        assert res == expect


def test_mapstr(configured_root):
    config = ConfigFile.from_yaml((configured_root / "config.yaml").as_posix())
    data = [
        ("with_cflavor", "bowtie2@custom_index:rRNA->STAR@default:genome", None),
    ]
    for name, mapstr, expect in data:
        mr, lr = mapstr_to_targets(mapstr)
        assert mr[0].input_name == "uBAM"
        assert mr[-1].input_name == "rRNA.bowtie2"
        assert lr[0].link_src == "genome.STAR"
        assert lr[0].link_name == "final"


def test_get_mapped_BAM_output(configured_root):
    config = ConfigFile.from_yaml((configured_root / "config.yaml").as_posix())
    project_df = ProjectDF(
        (configured_root / "project_df.csv").as_posix(), config=config
    )

    out_files = get_mapped_BAM_output(project_df=project_df, config=config)
    print(out_files)


def test_validation_cmdline_issue_54(configured_root):
    os.chdir(configured_root.as_posix())
    data = [
        ("flipped", "rRNA:bowtie2->genome:STAR", "test_hsa", True),
        ("species_missing", "bowtie2:rRNA->STAR:genome", "test_hs", False),
        (
            "with_cflavor",
            "rRNA:bowtie2@custom_index->genome:STAR@default",
            "test_hsa",
            True,
        ),
        (
            "unknown_cflavor",
            "rRNA:bowtie2@customX->genome:STAR@defaultBLA",
            "test_hsa",
            False,
        ),
        # ("flipped", "bowtie2:rRNA->STAR:genome", "rRNA:bowtie2->genome:STAR"),
    ]
    for name, mapstr, species, expect_pass in data:
        print(f"running test {name}")
        # add
        sm(
            "projects",
            "add-sample",
            "--project-id=test",
            f"--sample-id={name}",
            f"--map-strategy={mapstr}",
            f"--R1={spacemake_dir}/test_data/reads_chr22_R1.fastq.gz",
            f"--R2={spacemake_dir}/test_data/reads_chr22_R2.fastq.gz",
            f"--species={species}",
            expect_fail=not expect_pass,
        )
