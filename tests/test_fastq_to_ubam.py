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
        "--out-file",
        "/dev/null",
    )


def test_single():
    sm(
        "--read2",
        spacemake_dir + "test_data/reads_chr22_R2.fastq.gz",
        "--out-file",
        "/dev/null",
        """--cell='"A"'""",
    )


def test_preprocessing():
    pre = PreProcessor("quality:left=25,right=25")

    qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    sdata = SeqData(
        "@NS_mRNA:7 1:N:0:ACTGAGCG",
        "CCTGCTGGGAGGGG",
        "IIIIIIIIIIIIII",
        "CCTGCTGGGAGGGGGTGGGGGGAGGAGGAAGAGGTGGGGCTCTACTCTGATTAATTA",
        qual,
    )
    sdata = pre.process(sdata)
    assert sdata.q2 == qual
    assert len(sdata.tags) == 0

    qual = "###IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII###I#I#"
    sdata.q2 = qual
    res = pre.process(sdata)
    # print(res)
    tag_d = dict(res.tags)
    assert tag_d["T5"] == ["3"]
    assert tag_d["T3"] == ["7"]

    pre = PreProcessor("quality:left=25,right=25;polyA")

    sdata = SeqData(
        "@NS_mRNA:7 1:N:0:ACTGAGCG",
        "CCTGCTGGGAGGGG",
        "IIIIIIIIIIIIII",
        "CCTGCTGGGAGGGGGTGGGGGGAGGAGGAAGAGGTGGGGCTCTACTCTGATTAATTAAAAAAAAAGAGAAAAAAAAAAAAAAAAGGG",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIEEIEEIE#IIIIIIIIIIIIIIEIII###",
    )

    res = pre.process(sdata)
    tag_d = dict(res.tags)
    assert tag_d["T3"] == ["3", "28"]


if __name__ == "__main__":
    test_preprocessing()

# def test_minqual():
#     sm(
#         "--read2",
#         spacemake_dir + "test_data/reads_chr22_R2.fastq.gz",
#         "--out-bam",
#         "/dev/null",
#         "--min-qual",
#         "30",
#         """--cell='"A"'""",
#     )


def test_issue135():
    import spacemake.bin.fastq_to_uBAM as ubam
    from argparse import Namespace
    import io

    f1 = io.StringIO(
        "\n".join(
            [
                "@QNAME MUST NOT HAVE WHITESPACE",
                "CCTGCTGGGAGGGG",
                "+",
                "IIIIIIIIIIIIII",
                "",
            ]
        )
    )

    f2 = io.StringIO(
        "\n".join(
            [
                "@QNAME MUST NOT HAVE WHITESPACE",
                "CCTGCTGGGAGGGGGTGGGGGGAGGAGGAAGAGGTGGGGCTCTACTCTGATTAATTA",
                "+",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                "",
            ]
        )
    )

    pre = ubam.PreProcessor("quality:left=25,right=25;barcode")

    for sdata in ubam.SeqData.from_paired_end(f1, f2):
        sdata = pre.process(sdata)

        sam = sdata.render_SAM(flag=4)
        cols = sam.split()
        assert cols[0] == "QNAME"
        assert cols[1] == "4"
    # print(sam)


if __name__ == "__main__":
    test_issue135()
