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
    def build_from_barcodes(cls, src, l_prefix=10, l_suffix=15):
        # prefix index -> list index
        logger = logging.getLogger("BCIndex.build_from_barcodes")
        logger.debug("building temp set list")

        sets = [set() for i in range(2 ** (l_prefix * 2))]

        logger.debug("Ingesting sequences")
        n_seqs = 0
        rshift = 2 * l_suffix
        mask = np.uint64((1 << (2 * l_suffix)) - 1)
        T0 = time()

        for idx in src:
            n_seqs += 1
            # prefix = seq[:l_prefix]
            # suffix = seq[l_prefix:]

            # idx_p = _to_index(prefix)
            # # print(f"{prefix} -> {idx_p} L={2**(l_prefix*2)}")
            # idx_s = _to_index(suffix)

            idx_p = idx >> rshift
            idx_s = idx & mask

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
            l_prefix=l_prefix, l_suffix=l_suffix
        )  # just the last suffix we processed
        bci.PI = PI
        bci.SL = SL

        return bci

    def sanity_check(self):
        from collections import defaultdict

        logger = logging.getLogger("BCIndex.sanity_check")
        n_total = 0
        n_lists = 0
        nn_max = 0
        for idx_p in range(self.L_idx_p):
            l_ofs = self.PI[idx_p]
            if l_ofs > 0:
                n_lists += 1
                n = self.SL[l_ofs]  # number of suffixes to expect
                n_total += n
                nn_max = max(nn_max, n)

        logger.info(
            f"sanity check: {n_lists} prefix lists with total of {n_total} suffixes. max list size = {nn_max}"
        )

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


def make_insertions(idx, l=25):
    """
    Generate all single-base insertions/clipped to l again.
    Write into an array of l*4 uint32_t values.
    """

    full_mask = (1 << 2 * l)
    ext_mask = full_mask << 2
    full_mask -= 1
    ext_mask -= 1 # keep one extra base for insertions

    i = 0
    idx_l = idx << 2  # on base left-shifted

    print(f"{seqidx.uint64_to_seq(idx, l)} one base insertions")
    print(f"{seqidx.uint64_to_seq(idx_l, l)} left-shifted index")


    # print(f"{uint64_to_seq(idx, l)} one base insertions")
    # print(f"l={l} 1 << 2*l = {1 << (2*l):b}")
    # print(f"full_mask {full_mask:b}")

    i = 0
    idx_l = idx << 2  # one base left-shifted

    mask_l = ext_mask ^ 3 #<< (start * 2)) & ext_mask  # initial left mask
    mask_r = 0 #1 << (start * 2) - 1  # initial right mask
    # print(f"{uint64_to_seq(idx_l, l)} left-shifted index")
    for pos in range(l):
        # insertion at base pos
        idx0 = (idx_l & mask_l) | (idx & mask_r)
        
        # print(f"{mask_l:050b} mask_l")
        # print(f"{mask_r:050b} mask_r")
        # print(f"{uint64_to_seq(idx0, l)} IDX0 pos={pos} i={i}")

        for k in range(4):
            idx_var = idx0 | (k << pos*2)
            if pos < l - 1:
                print(f"{seqidx.uint64_to_seq(idx_var & full_mask, l)} pos={pos} k={k} i={i}")
                # variants[i] = idx_var & full_mask
                i += 1
            if pos > 0:
                print(f"{seqidx.uint64_to_seq(idx_var >> 2, l)} pos={pos} k={k} i={i}")
                # variants[i] = idx_var >> 2
                i += 1        
            # print(f"{uint64_to_seq(idx_var, l)} k={k} i={i}")

        mask_l = (mask_l << 2) & ext_mask
        mask_r = (mask_r << 2) | 3

    return i

    # for pos in range(l):
    #     # insertion at base pos
    #     mask_r = (1 << pos * 2) - 1
    #     mask_l = full_mask ^ ((mask_r << 2) | 11)
    #     idx0 = (idx_l & mask_l) | (idx_r & mask_r)

    #     for k in range(4):
    #         idx_var = idx0 | (k << pos * 2)
    #         i += 1

    #         print(f"{seqidx.uint64_to_seq(idx_var, l)} pos={pos} k={k} i={i}")



def make_deletions(idx, l=25):
    """
    Generate all single-base deletions, padded on either side to l again.
    Write into an array of l*4 uint32_t values.
    """

    variants = np.zeros(l * 8, dtype=np.uint64)
    full_mask = np.uint64(1)
    full_mask = full_mask << 2 * l
    full_mask -= 1

    print(f"{seqidx.uint64_to_seq(idx, l)} one base deletions")
    print(f"l={l} 1 << 2*l = {1 << (2*l):b}")
    print(f"full_mask {full_mask:b}")

    i = 0
    idx_l = idx >> 2  # last base deleted, shifted right

    mask_l = full_mask
    mask_r = 0
    print(f"{seqidx.uint64_to_seq(idx_l, l)} left-shifted index")
    for pos in range(l):
        # delete base pos
        idx0 = (idx_l & mask_l) | (idx & mask_r)

        # print(f"{mask_l:050b} mask_l")
        # print(f"{mask_r:050b} mask_r")
        print(f"{seqidx.uint64_to_seq(idx0, l-1)} IDX0 pos={pos} i={i}")

        for k in range(4):
            idx_var = (idx0 << 2) | k  # insert random base on the right
            variants[i] = idx_var
            i += 1
            print(f"{seqidx.uint64_to_seq(idx_var, l)} k={k} i={i} random right")

            idx_var = idx0 | (k << (2 * (l - 1)))  # insert random base on the left
            variants[i] = idx_var
            i += 1
            print(f"{seqidx.uint64_to_seq(idx_var, l)} k={k} i={i} random left")

        mask_l = (mask_l << 2) & full_mask
        mask_r = (mask_r << 2) | 3

    return i


def make_edit_dict(l=25):
    d = {
        0:'X', # no match
        1:'=', # exact match
    }
    
    i = 2
    
    # insertions first
    for pos in range(l):
        for b in "ACGT":
            if pos < l-1:
                d[i] = f"I{l-pos}{b}1"
                i += 1
            if pos > 0:
                d[i] = f"I{l-pos}{b}0"
                i += 1

    # substitutions
    for pos in range(l):
        for k in range(1, 4):
            d[i] = f"S{l-pos}+{k}"
            i += 1

    # deletions
    for pos in range(l):
        for b in "ACGT":
            d[i] = f"_{l-pos}..{b}"
            i += 1
            d[i] = f"_{b}..{l-pos}"
            i += 1
    
    return d


def testing():

    # freak = "AAAAAACAATATTAATGTGAGCTCG".encode("ascii")
    # f64 = seqidx.seq_to_uint64(freak)

    # n = make_insertions(f64)
    # print(f"generated {n} insertions")
    # #make_deletions(f64)
    # 1 / 0

    # bci = BCIndex.load_mmap(path=".", l_prefix=10, l_suffix=15)
    # # print(bci.query([freak]))

    # # print(bci.query_idx64([f64]))

    # freak = "AAAAAAACAATATTAATTGAGCTCG".encode("ascii")
    # f64 = seqidx.seq_to_uint64(freak)

    # ob1_hits = np.zeros(1, dtype=np.uint8)
    # seqidx.query_idx64_indel(
    #     [f64], ob1_hits, bci.PI, bci.SL, bci.l_prefix, bci.l_suffix
    # )
    # print(ob1_hits)
    # 1 / 0
    # # prefrix
    # bc = "ACGTACGTACGTACGTACGTACGTA"
    # idx = seqidx.seq_to_uint64(bytes(bc, "ascii"))
    # # shifts = np.array(50, dtype=np.uint64)
    # shifts = seqidx.make_shifts(idx, l=25)
    # print(f"original bc {bc} -> {idx}")
    # for s in shifts:
    #     sbc = seqidx.uint64_to_seq(s, 25)
    #     print(f"shifted versions {sbc} -> {s}")

    # 1 / 0
    # idx = _to_index(bc)
    # id2 = seqidx.seq_to_uint32(bc)
    # print(f"{bc} -> {idx} == {id2}")
    # # bc_ = _to_seq(idx, len(bc))
    # # print(f"{bc} -> {idx} -> {bc_}")
    # 1 / 0
    # test_data = list(reader("sample1k.txt.gz"))

    # test_data = sorted(set(reader("../longreads/fc_1_2_1414.txt.gz")))

    n_max = 100000000
    n_max = 1000000
    # n_max = 100000
    # n_max = 100
    # fname = "/data/rajewsky/home/zkliesm/ont_openst/reference/all_BCs/lib298_whitelist_allBCs.csv"
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

    # logging.debug("building BC index from test data")
    # T0 = time()
    # bci = BCIndex.build_from_barcodes(idx_data, l_prefix=10)
    # bci.sanity_check()
    # dT = time() - T0
    # rate = len(test_data) / dT / 1000
    # print(f"built BC index in {dT:.1f} seconds ({rate:.2f} k/sec)")

    logging.debug("loading stored index from mmap")
    bci = BCIndex.load_mmap(
        path=".", l_prefix=10, l_suffix=15
    )  # just the last suffix we processed


    logging.debug(
        f"testing d=1 neighbors of barcodes"
    )
    T0 = time()
    hit_variants = np.zeros(len(test_data), dtype=np.int16)
    hits = np.zeros(len(test_data), dtype=np.uint64)
    n_total_queries = seqidx.query_idx64_variants(
        list(test_data), hits, hit_variants, bci.PI, bci.SL, bci.l_prefix, bci.l_suffix
    )
    dT = time() - T0
    rate = len(test_data) / dT / 1000
    n_total_hits = (hit_variants > 0).sum()
    print(
        f"variants: {n_total_hits}/{len(test_data)} hits (match-rate = {n_total_hits/len(test_data):.4f}) in {dT:.1f} seconds ({rate:.2f} k/sec)"
    )
    print(
        f"total queries processed: {n_total_queries} ({n_total_queries/len(test_data):.2f} per barcode) rate: {n_total_queries / dT / 1000:.2f} k/sec)"
    )

    d = make_edit_dict(l=bci.l_prefix + bci.l_suffix)

    n_edits = np.bincount(hit_variants)

    for i, n in sorted(list(enumerate(n_edits)), key=lambda x: -x[1]):
        print(f"edit {i}: {d.get(i,'?')} -> {n} hits")

    for query, match, var_i in zip(test_data, hits, hit_variants):
        query_seq = seqidx.uint64_to_seq(query, bci.l_prefix + bci.l_suffix)
        match_seq = seqidx.uint64_to_seq(match, bci.l_prefix + bci.l_suffix)
        edit = d.get(var_i, "?")
        print(f"query {query_seq} matched {match_seq} via edit {edit}")

    # seqidx.query_idx64_indel(
    #     list(non_matched), ob1_hits, bci.PI, bci.SL, bci.l_prefix, bci.l_suffix
    # )
    # dT = time() - T0
    # rate = len(non_matched) / dT / 1000
    # n_uniq = (ob1_hits == 1).sum()
    # n_matches = (ob1_hits > 0).sum()
    # print(
    #     f"indel: {n_uniq}/{len(non_matched)} hits (uniq match-rate = {n_uniq/len(non_matched):.4f}) in {dT:.1f} seconds ({rate:.2f} k/sec)"
    # )
    # n_total_hits += n_uniq




    # bci.sanity_check()
    # logging.debug("testing completeness and correctness of stored barcode universe")
    # for seq, ref in zip(bci.dump(), test_data):
    #     if seq != ref:
    #         print(f"mismatch between {seq} and {ref}")

    # logging.debug("testing query with reference")
    # T0 = time()
    # hits = bci.query_idx64(test_data)
    # n_total_hits = hits.sum()
    # print(f"{hits.sum()} hits (match-rate = {n_total_hits/len(test_data):.4f})")
    # dT = time() - T0
    # rate = len(test_data) / dT / 1000
    # print(f"looked up {len(test_data)} in {dT:.1f} seconds ({rate:.2f} k/sec)")

    # # for idx in non_matched:
    # #     seq = seqidx.uint64_to_seq(idx, bci.l_prefix + bci.l_suffix)
    # #     print(f"non-matched barcode: {seq}")

    # # assert hits.all()

    # # non_matched = np.array(test_data)[~hits]
    # test_data = np.array(test_data)
    # non_matched = test_data[~hits]
    # if len(non_matched) > 0:
        # logging.debug("testing shifted versions of non-matched barcodes")
        # ob1_hits = np.zeros(len(non_matched), dtype=np.uint8)
        # seqidx.query_idx64_shifts(
        #     list(non_matched), ob1_hits, bci.PI, bci.SL, bci.l_prefix, bci.l_suffix
        # )
        # n_uniq = (ob1_hits == 1).sum()
        # n_matches = (ob1_hits > 0).sum()
        # print(
        #     f"shift: {n_uniq}/{n_matches} hits (uniq match-rate = {n_uniq/len(test_data):.4f})"
        # )
        # # bci.sanity_check()
        # n_total_hits += n_uniq
        # non_matched = non_matched[~ob1_hits] # still non-matched

        # lets create all off-by-one neighbors and see if they match
        # logging.debug(
        #     f"testing indel neighbors of {len(non_matched)} non-matched barcodes"
        # )
        # T0 = time()
        # ob1_hits = np.zeros(len(non_matched), dtype=np.uint8)
        # seqidx.query_idx64_indel(
        #     list(non_matched), ob1_hits, bci.PI, bci.SL, bci.l_prefix, bci.l_suffix
        # )
        # dT = time() - T0
        # rate = len(non_matched) / dT / 1000
        # n_uniq = (ob1_hits == 1).sum()
        # n_matches = (ob1_hits > 0).sum()
        # print(
        #     f"indel: {n_uniq}/{len(non_matched)} hits (uniq match-rate = {n_uniq/len(non_matched):.4f}) in {dT:.1f} seconds ({rate:.2f} k/sec)"
        # )
        # n_total_hits += n_uniq

        # logging.debug(
        #     f"testing off-by-one neighbors of {len(non_matched)} non-matched barcodes"
        # )
        # T0 = time()
        # ob1_hits = np.zeros(len(non_matched), dtype=np.uint8)
        # seqidx.query_idx64_off_by_one(
        #     list(non_matched), ob1_hits, bci.PI, bci.SL, bci.l_prefix, bci.l_suffix
        # )
        # dT = time() - T0
        # rate = len(non_matched) / dT / 1000
        # n_uniq = (ob1_hits == 1).sum()
        # n_matches = (ob1_hits > 0).sum()
        # print(
        #     f"hamming1: {n_uniq}/{len(non_matched)} hits (uniq match-rate = {n_uniq/len(non_matched):.4f}) in {dT:.1f} seconds ({rate:.2f} k/sec)"
        # )
        # n_total_hits += n_uniq

    print(
        f"total matched barcodes: {n_total_hits}/{len(test_data)} ({n_total_hits/len(test_data):.4f})"
    )


if __name__ == "__main__":
    testing()
