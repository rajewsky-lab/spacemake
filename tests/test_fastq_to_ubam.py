import multiprocessing as mp

import pytest
import sys
import os
from spacemake.bin.fastq_to_uBAM import *


spacemake_dir = os.path.dirname(__file__) + "/../"


@pytest.fixture(scope="session")
def test_root(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("root")
    sm_path = os.path.dirname(__file__)
    # make a tmp-copy of the test_config.yaml
    # def_config = os.path.join(sm_path, "../test_data/test_config.yaml")
    # os.system(f"cp {def_config} {tmp / 'config.yaml'}")

    # test_pdf =  os.path.join(sm_path, "../test_data/test_project_df.csv")
    # os.system(f"cp {test_pdf} {tmp / 'project_df.csv'}")

    return tmp


def sm(*argc, expect_fail=False):
    sys.argv = [
        "fastq_to_uBAM.py",
    ] + list(argc)
    res = cmdline()
    print("got result", res)
    from spacemake.errors import SpacemakeError

    if expect_fail:
        assert isinstance(res, SpacemakeError) == True
    else:
        assert isinstance(res, Exception) == False


def test_help():
    try:
        sm("--help")
    except SystemExit:
        pass


def test_dropseq():
    sm(
        "--read1",
        spacemake_dir + "test_data/reads_chr22_R1.fastq.gz",
        "--read2",
        spacemake_dir + "test_data/reads_chr22_R2.fastq.gz",
        "--out-bam",
        "/dev/null",
    )


def test_single():
    sm(
        "--read2",
        spacemake_dir + "test_data/reads_chr22_R2.fastq.gz",
        "--out-bam",
        "/dev/null",
        """--cell='"A"'""",
    )


def test_minqual():
    sm(
        "--read2",
        spacemake_dir + "test_data/reads_chr22_R2.fastq.gz",
        "--out-bam",
        "/dev/null",
        "--min-qual",
        "30",
        """--cell='"A"'""",
    )


def test_issue135():
    import spacemake.bin.fastq_to_uBAM as ubam
    from argparse import Namespace

    args = Namespace(
        bam_tags="CR:{cell},CB:{cell},MI:{UMI},RG:A",
        min_len=18,
        min_qual_trim=0,
        cell="r1[8:20][::-1]",
        UMI="r1[0:8]",
        seq="r2",
        qual="r2_qual",
        disable_safety=False,
    )

    fmt = ubam.make_formatter_from_args(args)
    attrs = fmt(
        r2_qname="QNAME MUST NOT HAVE WHITESPACE",
        r1="ACGTACGT",
        r1_qual="########",
        r2="TGCATGCATGCATGCA",
        r2_qual="################",
    )
    sam = ubam.make_sam_record(flag=4, **attrs)
    cols = sam.split()
    assert cols[0] == "QNAME"
    assert cols[1] == "4"
    # print(sam)


# if __name__ == "__main__":
#     test_issue135()
