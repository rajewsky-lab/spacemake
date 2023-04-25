import pytest
import os
from spacemake.annotation import GenomeAnnotation

@pytest.fixture(scope="session")
def test_compile(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("simple")
    spacemake_path = os.path.dirname(__file__)
    print("COMPILING ANNOTATION")
    ga = GenomeAnnotation.from_GTF(
        os.path.join(spacemake_path, "../test_data/simple_annotation.gtf"),
        df_cache=(tmp / "simple.tsv").as_posix()
    )
    _ = ga.compile((tmp / "compiled/").as_posix())
    return tmp / "compiled"

def test_single_gene(test_compile):
    regions = [
        (('chr1', '+', [(0, 74), ]), ("A", "N", "c")),
        (('chr1', '+', [(75, 150), ]), ("A", "C", "c")),
        (('chr1', '+', [(0, 150), ]), ("A", "N,C", "c")),
        (('chr1', '+', [(75, 151), ]), ("A", "C,I", "c")),
        (('chr1', '-', [(0, 150), ]), ("A", "n", "c")),
    ]
    print(test_compile)
    ga = GenomeAnnotation.from_compiled_index(test_compile)
    
    for (chrom, strand, blocks), expect in regions:
        gn, gf, gt = ga.get_annotation_tags(chrom, strand, blocks)
        print(f'gn "{gn}"; gf "{gf}"; gt "{gt}"; expect={expect}')


