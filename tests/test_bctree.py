import pytest
import os
import spacemake.bctree as bt

@pytest.fixture(scope="session")
def test_compile(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("bctree")
    spacemake_path = os.path.dirname(__file__)
    compiled_path = tmp / "bctree.npy"

    bc = bt.compile([os.path.join(spacemake_path, "bc_test_data_0001.txt.gz")], dbname=compiled_path.as_posix())
    return compiled_path

@pytest.fixture(scope="session")
def test_load(test_compile):
    bc = bt.BCTree.from_disk(test_compile.as_posix())
    return bc

def test_batch(test_load):
    bc = test_load
    bc.check_contains_batch()

def test_hull(test_load):
    bc = test_load
    
    # edit distance 1 (aka "hull") generation benchmark
    for idx in bc.buf['idx'][:bc.n]:
        bt.bctree.hull(idx, 24)
    





