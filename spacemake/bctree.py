import numpy as np
import numpy
import pyximport
pyximport.install(setup_args={"include_dirs": numpy.get_include()}, reload_support=True)
import spacemake.cython.bctree as bctree
from time import time
import logging

LEAF = 0xFFFFFFFE

class BCTree:
    logger = logging.getLogger("spacemake.BCTree")
    dtype = [('idx', np.uint64), ('left', np.uint32), ('right', np.uint32), ('x', np.uint16), ('y', np.uint16), ('tile', np.uint16)]

    def __init__(self, N=int(2E9)):
        self.N = N
        self.n = 0
        self.buf = np.zeros(N, dtype=self.dtype)
        self.root = None

    def load_tile_data(self, fnames):
        import gzip
        t0 = time()
        n = 0
        for fname in fnames:
            tile = int(fname.split(".txt")[0].split('_')[-1])
            if fname.endswith(".gz"):
                src = gzip.open(fname, 'rt')
            else:
                src = open(fname)

            n_read = bctree.read_txt_to_buf(src, self.buf[n:], tile=tile)
            n += n_read
        
        self.n = n
        t1 = time()
        self.logger.info(f"{1000 * (t1-t0):.2f} msec loading and converting {self.n} barcode records")
        
        # truncate array
        self.N = self.n
        self.buf = self.buf[:n]
    
    def sort(self):
        # TODO: replace with faster, parallel implementation (e.g. parallel numpy)
        try:
            import pnumpy as pn
            pn.thread_setworkers(16)
        except ImportError:
            self.logger.warning("Can not import pnumpy (parallel numpy). Sorting will be slow.")

        t0 = time()
        
        ## in-place sorting is slow
        #self.buf[:self.n].sort()
        
        # faster to sort on the idx integers only and then re-order
        I = self.buf['idx'][:self.n].argsort()
        self.buf = self.buf[I]

        t1 = time()
        dt = t1 - t0
        self.logger.info(f"{1000 * dt:.2f} msec sorting {self.n} barcode records")

    def build_BST(self):
        t0 = time()
        
        # build-tree
        self.root = bctree.build_tree(self.buf, 0, self.n)
        
        t1 = time()
        dt = t1 - t0
        self.logger.info(f"{1000 * dt:.2f} msec building {self.n} barcode BST")

    def check_contains_batch(self):
        assert self.root is not None
        t0 = time()

        # checking in batch
        bctree.check_contains_batch(self.buf, self.root, self.buf['idx'][:self.n])

        t1 = time()
        dt = t1 - t0
        self.logger.info(f"{1000 * dt:.2f} msec batch checking for match of {self.n} existing barcodes ({0.001 * self.n/dt:.2f} kBC/sec)")
        
    def save_to_disk(self, filename):
        np.save(filename, self.buf)

    @classmethod
    def from_disk(cls, filename):
        buf = np.lib.format.open_memmap(filename, dtype=cls.dtype, mode='r')
        bc = cls(N=len(buf))
        bc.n = len(buf)
        bc.buf = buf
        bc.root = bc.n >> 1

        bc.logger.info(f"loaded BCTree data for {bc.n} barcodes from file '{filename}'.")
        return bc

def compile(fnames, dbname='bctree.npy'):
    bc = BCTree()
    bc.load_tile_data(fnames)
    bc.sort()
    bc.build_BST()

    bc.save_to_disk(dbname)
    return bc


#if __name__ == "__main__":
    #logging.basicConfig(level=logging.INFO)

    ## 2.4M smallest example. Should be very fast
    # bc = compile(["fc_010_L3_tile_2267.txt"])

    ## 8M barcodes
    # bc = compile(["fc_010_L3_tile_2267.txt", "fc_010_L3_tile_2268.txt", "fc_010_L3_tile_2269.txt"], 'bctree.npy')
    
    ## 216M barcodes + testing gz input -> takes ~4.5 GB on disk
    # from glob import glob
    # fnames = glob("/data/rajewsky/projects/slide_seq/puck_data/seq_scope_inhouse_fc_sts_nova5_FC_010_230109/fc_010_L3_tile_22*.gz")
    # print(fnames)
    # bc = compile(fnames)

    # bc.check_contains_batch()
    
    # bc = BCTree.from_disk('bctree.npy')
    # bc.check_contains_batch()

    # t5 = time()
    # # edit distance 1 (aka "hull") generation benchmark
    # for idx in bc.buf['idx'][:bc.n]:
    #     bctree.hull(idx, 24)

    # t6 = time()
    # dt = t6 - t5
    # print(f"{1000 * dt:.2f} msec hull generation for {bc.n} existing barcodes ({0.001 * bc.n/dt:.2f} kBC/sec)")


    # t5 = time()
    # # edit distance 1 (aka "hull") generation benchmark
    # for idx in bc.buf['idx'][:bc.n]:
    #     bctree.hull(idx, 24)

    # t6 = time()
    # dt = t6 - t5
    # print(f"{1000 * dt:.2f} msec hull generation for {bc.n} existing barcodes ({0.001 * bc.n/dt:.2f} kBC/sec)")
