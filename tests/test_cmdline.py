import pytest
import sys
import os
from spacemake.cmdline import *

from fixtures import initialized_root, with_species, configured_root, sm, spacemake_dir


def test_init(initialized_root):
    assert os.path.exists(initialized_root.as_posix() + '/config.yaml')


def test_parsers(configured_root):
    get_project_sample_parser()
    get_add_sample_sheet_parser()
    
    parser = get_sample_main_variables_parser()
    get_action_sample_parser(parser.add_subparsers(), "add", lambda *argc, **kw: None)

    get_sample_extra_info_parser()
    get_data_parser()
    get_run_parser()

    make_main_parser()


def test_subcommands(configured_root):
    # test everything after init
    os.chdir(configured_root.as_posix())

    sm("projects", "list")
    sm("config", "list_adapter_flavors")


def test_config_barcodes(initialized_root):
    os.chdir(initialized_root.as_posix())

    # add
    sm(
        "config", "add_barcode_flavor",
        "--name", "fc_sts_miniseq",
        "--umi", "r2[0:9]", "--cell_barcode", "r1[2:27]"
    )
    # add same bc flavor twice? -> error
    sm(
        "config", "add_barcode_flavor",
        "--name", "fc_sts_miniseq",
        "--umi", "r2[0:9]", "--cell_barcode", "r1[2:27]",
        expect_fail=True
    )
    # delete
    sm(
        "config", "delete_barcode_flavor",
        "--name", "fc_sts_miniseq",
    )
    # re-add
    sm(
        "config", "add_barcode_flavor",
        "--name", "fc_sts_miniseq",
        "--umi", "r2[0:8]", "--cell_barcode", "r1[2:27]"
    )
    # update
    sm(
        "config", "update_barcode_flavor",
        "--name", "fc_sts_miniseq",
        "--umi", "r2[0:9]"
    )


def test_config_adapters(initialized_root):
    os.chdir(initialized_root.as_posix())
    # add
    sm(
        "config", "add_adapter",
        "--name", "testy", "--seq", "ACGTACGTACGTACGT",
    )
    sm(
        "config", "add_adapter",
        "--name", "testy", "--seq", "ACGTACGTACGTACGT",
        expect_fail=True
    )
    # update
    sm(
        "config", "update_adapter",
        "--name", "testy", "--seq", "ACGTACGTACGTACGTA",
    )
    # delete
    sm(
        "config", "delete_adapter",
        "--name", "testy"
    )
    # re-add
    sm(
        "config", "add_adapter",
        "--name", "testy", "--seq", "ACGTACGTACGTACGTACGT",
    )


def test_config_adapter_flavors(initialized_root):
    os.chdir(initialized_root.as_posix())

    test_flavor = (
        "--name", "testy1",
        "--cut_left", 
        "SMART:min_overlap=10:max_errors=0.1",
        "--cut_right", 
        "Q:min_base_qual=30", 
        "polyG:min_overlap=3:max_errors=0.2",
        "polyA:min_overlap=3:max_errors=0.25"
    )

    # add
    sm("config", "add_adapter_flavor", *test_flavor)
    sm("config", "add_adapter_flavor", *test_flavor, expect_fail=True)
    # update -> not implemented
    sm("config", "update_adapter_flavor", "--name", "testy1", expect_fail=True)
    # delete
    sm("config", "delete_adapter_flavor", "--name", "testy1")
    # re-add
    sm("config", "add_adapter_flavor", *test_flavor)


def test_species(with_species):
    pass


def test_puck(initialized_root):
    os.chdir(initialized_root.as_posix())
    test_puck = (
        "--name=test_puck",
        "--width_um=2",
        "--spot_diameter_um=0.5",
    )
    # add
    sm("config", "add_puck", *test_puck)
    sm("config", "add_puck", *test_puck, expect_fail=True)
    # update
    sm("config", "update_puck", "--name=test_puck", "--width_um=3")
    # delete
    sm("config", "delete_puck", "--name=test_puck")
    # re-add
    sm("config", "add_puck", *test_puck)


def test_runmode(initialized_root):
    os.chdir(initialized_root.as_posix())
    # add
    sm("config", "add_run_mode", "--name=spatial_rm", "--umi_cutoff=10")
    sm("config", "add_run_mode", "--name=spatial_rm2","--umi_cutoff=10")
    # delete
    sm("config", "delete_run_mode", "--name=spatial_rm2")
    # edit
    sm("config", "update_run_mode", "--name=spatial_rm", "--umi_cutoff=1")
    from spacemake.config import get_global_config
    config = get_global_config()
    config.dump()

def test_sample(configured_root):
    os.chdir(configured_root.as_posix())
    test_sample = (
        f'--R1={spacemake_dir}/test_data/reads_chr22_R1.fastq.gz', 
        f'--R2={spacemake_dir}/test_data/reads_chr22_R1.fastq.gz', 
        '--map_strategy=genome:STAR:final',
        '--species=test_hsa'
    )

    # add
    sm(
        "projects", "add_sample",
        "--project_id=test", "--sample_id=test1", 
        *test_sample
    )
    # delete
    sm(
        "projects", "delete_sample",
        "--project_id=test", "--sample_id=test1", 
    )
    # re-add
    sm(
        "projects", "add_sample",
        "--project_id=test", "--sample_id=test1", 
        *test_sample
    )
    # update
    sm(
        "projects", "update_sample",
        "--project_id=test", "--sample_id=test1",
        '--map_strategy=rRNA:bowtie2->genome:STAR'
    )

def test_fill_project_df(with_species):
    os.chdir(with_species.as_posix())
    test_project_data = [
        (
            "test_hsa",
            "test",
            "test_01",
            f"{spacemake_dir}/test_data/reads_chr22_R1.fastq.gz",
            f"{spacemake_dir}/test_data/reads_chr22_R2.fastq.gz",
            "--map_strategy=genome:STAR:final",
        ),
        (
            "test_hsa",
            "test",
            "test_02.2",
            f"{spacemake_dir}/test_data/reads_chr22_R1.fastq.gz",
            f"{spacemake_dir}/test_data/reads_chr22_R2.fastq.gz",
            "--map_strategy=rRNA:bowtie2->miRNA:bowtie2->genome:STAR:final",
        ),
        (
            "test_hsa",
            "test",
            "test_03_nofinal",
            f"{spacemake_dir}/test_data/reads_chr22_R1.fastq.gz",
            f"{spacemake_dir}/test_data/reads_chr22_R2.fastq.gz",
            "--map_strategy=rRNA:bowtie2->miRNA:bowtie2->genome:STAR",
        ),
        (
            "test_hsa",
            "test",
            "test_bulk",
            "None",
            f"{spacemake_dir}/test_data/reads_chr22_R2.fastq.gz",
            "--map_strategy=rRNA:bowtie2->miRNA:bowtie2->genome:STAR:final"
            " --barcode_flavor=visium",
        ),
        (
            "test_hsa",
            "tile",
            "tile_1",
            f"{spacemake_dir}/test_data/reads_chr22_R1.fastq.gz",
            f"{spacemake_dir}/test_data/reads_chr22_R2.fastq.gz",
            (
                "--map_strategy=rRNA:bowtie2->miRNA:bowtie2->genome:STAR:final"
                f" --puck_barcode_file {spacemake_dir}/test_data/tile_1.txt"
                " --puck slide_seq --run_mode slide_seq"
            )
        ),
        (
            "test_hsa",
            "tile",
            "tile_2",
            f"{spacemake_dir}/test_data/reads_chr22_R1.fastq.gz",
            f"{spacemake_dir}/test_data/reads_chr22_R2.fastq.gz",
            (
                "--map_strategy=rRNA:bowtie2->miRNA:bowtie2->genome:STAR:final"
                f" --puck_barcode_file {spacemake_dir}/test_data/tile_2.txt"
                " --puck slide_seq --run_mode slide_seq"
            )
        ),
        (
            "test_hsa",
            "tile",
            "tile_both",
            f"{spacemake_dir}/test_data/reads_chr22_R1.fastq.gz",
            f"{spacemake_dir}/test_data/reads_chr22_R2.fastq.gz",
            (
                "--map_strategy=rRNA:bowtie2->miRNA:bowtie2->genome:STAR:final"
                f" --puck_barcode_file {spacemake_dir}/test_data/tile_1.txt {spacemake_dir}/test_data/tile_2.txt"
                " --puck slide_seq --run_mode slide_seq"
            )
        ),
    ]
    for species, project_id, sample_id, R1, R2, extra in test_project_data:
        # add
        sm(
            "projects", "add_sample",
            f"--project_id={project_id}", f"--sample_id={sample_id}",
            f"--species={species}",
            f"--R1={R1}", f"--R2={R2}", *extra.split(' ')
        )


def test_run(configured_root):
    # test everything after init
    os.chdir(configured_root.as_posix())
    sm("run", "-np", "--cores=8")