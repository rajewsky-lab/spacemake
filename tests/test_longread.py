import pytest
import os
from fixtures import initialized_root, with_species, sm, spacemake_dir

def test_mm2_junc_bed(with_species):
    """Test the --junc-bed wiring to mm2 mapping"""
    os.chdir(with_species.as_posix())
    sm(
        "projects",
        "add-sample",
        "--project-id=test",
        "--sample-id=test_mm2_junc_bed",
        f"--R1={spacemake_dir}/test_data/simple.reads1.fastq.gz",
        f"--R2={spacemake_dir}/test_data/simple.reads2.fastq.gz",
        "--map-strategy=mm2:genome:final",
        "--species=test_hsa",
    )
    sm("run", "-np", "--cores=8")
