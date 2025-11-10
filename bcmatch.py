import re
import sys
import numpy as np
import seqidx
from seqidx import seq_to_uint32 as _to_index

# def _to_index(seq):
#     assert len(seq) <= 16  # else uint32 has not enough bits
#     bits = {
#         "a": np.uint32(0),
#         "c": np.uint32(1),
#         "g": np.uint32(2),
#         "t": np.uint32(3),
#         "A": np.uint32(0),
#         "C": np.uint32(1),
#         "G": np.uint32(2),
#         "T": np.uint32(3),
#     }

#     idx = np.uint32(0)
#     i = 0
#     for nt in seq[::-1]:
#         idx |= bits[nt] << i
#         i += 2

#     return idx


def _to_seq(idx, l):
    bases = {0: "A", 1: "C", 2: "G", 3: "T"}
    seq = [""] * l
    for i in range(l):
        seq[l - i - 1] = bases[idx & 3]
        idx = idx >> 2

    return "".join(seq)


import logging
from time import time

logging.basicConfig(level=logging.DEBUG)


class BCIndex:
    def __init__(self, l_prefix=10, l_suffix=15):
        self.l_prefix = l_prefix
        self.l_suffix = l_suffix
        self.L_idx_p = 2 ** (l_prefix * 2)

    @classmethod
    def load_mmap(cls, path=".", l_prefix=10, l_suffix=15, **kw):
        bci = cls(l_prefix=l_prefix, l_suffix=l_suffix)
        bci.PI = np.memmap(
            f"{path}/PI.mmap", mode="r", shape=(2 ** (l_prefix * 2)), dtype=np.uint32
        )
        bci.SL = np.memmap(
            f"{path}/SL.mmap", mode="r", shape=2**30, dtype=np.uint32
        )  # 1E9
        return bci

    @classmethod
    def build_from_barcodes(cls, src, l_prefix=10, **kw):
        # prefix index -> list index
        logger = logging.getLogger("BCIndex.build_from_barcodes")
        logger.debug("building temp set list")

        sets = [set() for i in range(2 ** (l_prefix * 2))]

        logger.debug("Ingesting sequences")
        n_seqs = 0
        T0 = time()
        for seq in src:
            n_seqs += 1
            prefix = seq[:l_prefix]
            suffix = seq[l_prefix:]

            idx_p = _to_index(prefix)
            # print(f"{prefix} -> {idx_p} L={2**(l_prefix*2)}")
            idx_s = _to_index(suffix)

            sets[idx_p].add(idx_s)

            if n_seqs % 1000000 == 0:
                dT = time() - T0
                rate = n_seqs / dT / 1000
                logger.info(
                    f"processed {n_seqs} sequences in {dT:.1f} seconds ({rate:.1f}k/sec)"
                )

        dT = time() - T0
        rate = n_seqs / dT / 1000
        logger.info(
            f"Finished ingest. Processed {n_seqs} sequences in {dT:.1f} seconds ({rate:.1f}k/sec)"
        )

        logger.debug("allocating buffers")
        # prefix index -> suffix-list starting offset
        PI = np.memmap(
            "PI.mmap", mode="w+", shape=(2 ** (l_prefix * 2)), dtype=np.uint32
        )
        # space for suffix-lists
        SL = np.memmap("SL.mmap", mode="w+", shape=2**30, dtype=np.uint32)  # 1E9

        logger.debug("packing suffix lists")
        ofs = 1  # one buffer byte, and keep 0 as "no list"
        for idx_p in range(2 ** (l_prefix * 2)):
            S = sets[idx_p]
            if not S:
                continue

            # compile the suffix set into a packed suffix list
            sl = sorted(sets[idx_p])
            L = len(sl)
            PI[idx_p] = ofs
            SL[ofs] = L
            ofs += 1
            SL[ofs : ofs + L] = sl  # assign values
            ofs += L

        logger.debug(f"needed {ofs *4} bytes for packed suffix lists")

        bci = cls(
            l_prefix=l_prefix, l_suffix=len(suffix)
        )  # just the last suffix we processed
        bci.PI = PI
        bci.SL = SL

        return bci

    def dump(self):
        for idx_p in np.arange(self.L_idx_p):
            l_ofs = self.PI[idx_p]
            if l_ofs > 0:
                prefix = _to_seq(idx_p, self.l_prefix)
                n = self.SL[l_ofs]  # number of suffixes to expect
                l_ofs += 1
                for i in range(n):
                    idx_s = self.SL[l_ofs + i]
                    suffix = _to_seq(idx_s, self.l_suffix)
                    yield f"{prefix}{suffix}"

    def query(self, bc_list):
        logger = logging.getLogger("BCIndex.query()")
        T0 = time()
        hits = np.zeros(len(bc_list), dtype=np.bool)
        seqidx.query(bc_list, hits, self.PI, self.SL, self.l_prefix)
        # for i, bc in enumerate(bc_list):
        #     prefix = bc[: self.l_prefix]
        #     idx_p = _to_index(prefix)
        #     ofs = self.PI[idx_p]
        #     # logging.debug(f"testing {bc}: {prefix} -> {idx_p} -> {ofs}")
        #     if ofs > 0:
        #         suffix = bc[self.l_prefix :]
        #         idx_s = _to_index(suffix)
        #         # look for suffix in suffix list
        #         n = self.SL[ofs]
        #         ofs += 1

        #         # TODO: replace with bisect for lookup time from O(N) to O(log(N))
        #         for j in range(n):
        #             # logging.debug(f"{suffix} -> {idx_s} == {self.SL[ofs]}?")
        #             if self.SL[ofs] == idx_s:
        #                 hits[i] = True
        #                 # logging.debug(f"hit found")
        #                 break
        #             ofs += 1

        dT = time() - T0
        rate = len(bc_list) / dT / 1000
        logger.debug(f"queried {len(bc_list)} in {dT:.1f} seconds ({rate:.2f} k/sec)")

        return hits

    def query_idx64(self, bc_list):
        logger = logging.getLogger("BCIndex.query_idx64()")
        T0 = time()
        hits = np.zeros(len(bc_list), dtype=np.bool)
        seqidx.query_idx64(
            bc_list, hits, self.PI, self.SL, self.l_prefix, self.l_suffix
        )
        dT = time() - T0
        rate = len(bc_list) / dT / 1000
        logger.debug(f"queried {len(bc_list)} in {dT:.1f} seconds ({rate:.2f} k/sec)")

        return hits


def reader(fname, n_max=None):
    from time import time

    n = 0
    if fname.endswith(".gz"):
        logging.debug("opening gzip file")
        import isal.igzip_threaded as igzip

        _open = igzip.open
    else:
        _open = open

    T0 = time()
    for line in _open(fname, "rt"):
        if line.startswith("cell_bc"):
            continue
        n += 1
        if n_max and n > n_max:
            break

        seq = line.rstrip().split("\t")[0].replace("N", "A")
        yield bytes(seq, "ascii")

    dT = time() - T0
    rate = n / dT / 1000
    logging.debug(f"read {n} barcodes in {dT:.1f} seconds ({rate:.2f} k/sec)")


def testing():

    # # prefrix
    # bc = "ACGTACGTAC"
    # idx = _to_index(bc)
    # id2 = seqidx.seq_to_uint32(bc)
    # print(f"{bc} -> {idx} == {id2}")
    # # bc_ = _to_seq(idx, len(bc))
    # # print(f"{bc} -> {idx} -> {bc_}")
    # 1 / 0
    # test_data = list(reader("sample1k.txt.gz"))

    # test_data = sorted(set(reader("../longreads/fc_1_2_1414.txt.gz")))

    n_max = 100000000
    # fname = "/data/local/rajewsky/home/zkliesm/ont_openst/reference/all_BCs/lib298_whitelist_allBCs.csv"
    fname = "bc_to_match.txt"
    logging.debug("loading test data as uint64")
    idx_data = seqidx.load_and_unique_sorted_barcodes(
        fname,
        n_max=n_max,
        unique=False,
        # buf_size=1024 * 16,
    )
    # logging.debug("loading test data as strings")
    # str_data = list(reader(fname, n_max=n_max))
    # logging.debug("comparing string and uint64 representations")
    # print(f"strings {len(str_data)} indices {len(idx_data)}")
    # for i, (sbc, idx) in enumerate(zip(str_data, idx_data)):
    #     _sbc = seqidx.uint64_to_seq(idx, len(sbc))
    #     if sbc.decode("ascii") != _sbc:
    #         print(f"{i} : mismatch between {sbc} and {_sbc}")

    test_data = idx_data

    # for i in range(10):
    #     idx = test_data[i]
    #     bc = seqidx.uint64_to_seq(idx, 25)

    #     prefix = bc[:10]
    #     suffix = bc[10:]

    #     _pre = seqidx.uint64_to_seq(idx >> (2 * 15), 10)
    #     _suf = seqidx.uint64_to_seq(idx & ((1 << (2 * 15)) - 1), 15)

    #     logging.debug(
    #         f"test data sample {i}: {idx} {bc} prefix={prefix} suffix={suffix} ->_pre={_pre} _suf={_suf}"
    #     )
    # dT = time() - T0
    # rate = len(test_data) / dT / 1000
    # print(f"looked up {len(test_data)} in {dT:.1f} seconds ({rate:.2f} k/sec)")
    # assert hits.all()

    # 1 / 0
    logging.debug("loading stored index from mmap")
    bci = BCIndex.load_mmap(
        path=".", l_prefix=10, l_suffix=15
    )  # just the last suffix we processed
    # logging.debug("testing completeness and correctness of stored barcode universe")
    # for seq, ref in zip(bci.dump(), test_data):
    #     if seq != ref:
    #         print(f"mismatch between {seq} and {ref}")

    logging.debug("testing query with reference")
    hits = bci.query_idx64(test_data)

    print(f"{hits.sum()} hits (match-rate = {hits.sum()/len(test_data):.4f})")
    # assert hits.all()


if __name__ == "__main__":
    testing()
