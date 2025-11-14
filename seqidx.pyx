# cython: boundscheck=False, wraparound=False, cdivision=True, nonecheck=False, initializedcheck=False
## cython: language_level=3
import cython
cimport cython
import numpy as np
# from libc.stdint cimport uint8_t
from libc.stdint cimport uint64_t, uint32_t, uint8_t
from numpy cimport ndarray, uint32_t as np_uint32_t, uint8_t as np_uint8_t
from cython cimport Py_ssize_t
from cpython.bytes cimport PyBytes_AsStringAndSize

from types import *
from libc cimport stdlib, stdio
from libc.string cimport strncmp


cdef uint32_t[256] base_map
for i in range(256):
    base_map[i] = 0

# unknown bases map to 0 (A)
base_map[ord('C')] = 1
base_map[ord('G')] = 2
base_map[ord('T')] = 3
base_map[ord('c')] = 1
base_map[ord('g')] = 2
base_map[ord('t')] = 3

cdef uint64_t[16] dimers
for i in range(16):
    dimers[i] = i

cpdef uint32_t seq_to_uint32(bytes seq):
    """
    Encode up to 16 bases (A,C,G,T) into a 32-bit unsigned integer.
    Packing is left-to-right: first base occupies the highest-used bits.
    """
    cdef Py_ssize_t n = len(seq)
    if n > 16:
        raise ValueError("Sequence too long (max 16 bases)")

    cdef uint32_t result = 0
    cdef unsigned char nt
    
    for nt in seq:
        result <<= 2
        result |= base_map[nt]

    return result

@cython.cfunc
@cython.inline
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
cdef uint32_t seq_to_uint32_from_buf_nogil(const unsigned char *buf, Py_ssize_t start, Py_ssize_t length) nogil:
    """Compute uint32 code for a raw bytes buffer without any Python indexing (nogil)."""
    cdef Py_ssize_t k, end = start + length
    cdef uint32_t result = 0
    for k in range(start, end):
        result <<= 2
        result |= base_map[buf[k]]

    return result

cpdef uint64_t seq_to_uint64(bytes seq):
    """
    Encode up to 32 bases (A,C,G,T) into a 64-bit unsigned integer.
    Packing is left-to-right: first base occupies the highest-used bits.
    """
    cdef Py_ssize_t n = len(seq)
    if n > 32:
        raise ValueError("Sequence too long (max 32 bases)")

    cdef uint64_t result = 0
    cdef unsigned char nt
    
    for nt in seq:
        result <<= 2
        result |= base_map[nt]

    return result

@cython.cfunc
@cython.inline
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
cdef uint64_t seq_to_uint64_from_buf_nogil(const unsigned char *buf, Py_ssize_t start, Py_ssize_t length) nogil:
    """Compute uint64 code for a raw bytes buffer without any Python indexing (nogil)."""
    # assert length <= 32
    cdef Py_ssize_t k, end = start + length
    cdef uint64_t result = 0
    for k in range(start, end):
        result <<= 2
        result |= base_map[buf[k]]

    return result


cpdef str uint64_to_seq(uint64_t idx, int k):
    """
    Decode a 64-bit unsigned integer into a DNA sequence of length k (max 32).
    """
    cdef char[32] buf
    cdef char[4] bases = ['A', 'C', 'G', 'T']
    cdef int i
    cdef uint64_t bits

    if k > 32:
        raise ValueError("k too large (max 32)")

    for i in range(k):
        buf[k - i - 1] = bases[idx & 0b11]
        idx >>= 2

    return bytes(buf[:k]).decode('ascii')




def load_and_unique_sorted_barcodes(fname, int k=25, int n_max=0, bint unique=False, size_t buf_size=4096):
    """
    Load barcodes from a text file (one per line, tab-separated fields).
    If n_max > 0, stop after reading n_max barcodes.
    """
    from time import time
    import logging

    cdef int n = 0
    cdef ssize_t n_read
    cdef size_t _bs 
    cdef uint64_t idx64

    cdef set uniq = set()
    cdef list bc_list = []
    T0 = time()

    # bypass python file-io for low-level C
    cdef char* buffer = <char*>stdlib.malloc(4*buf_size) # file I/O buffer
    cdef char* line = <char*>stdlib.malloc(buf_size) # max line length 4k
    cdef stdio.FILE *fin = stdio.fopen(fname.encode('utf-8'), 'r')
    stdio.setvbuf(fin, buffer, stdio._IOFBF, buf_size)

    while(True):
        _bs = buf_size
        n_read = stdio.getline(&line, &_bs, fin)
        if n_read <= 0:
            break
        
        if n_read < k:
            # skip too short barcodes
            continue
    #     if line.startswith("cell_bc"):
    #         continue

        # if (n < 2) and (strncmp(buffer, b"cell_bc", 7) == 0):
        #     continue

        n += 1
        if n_max and n > n_max:
            break

        idx64 = seq_to_uint64_from_buf_nogil(<const unsigned char*>line, 0, k)

        if unique:
            if not idx64 in uniq:
                uniq.add(idx64)
                bc_list.append(idx64)
        else:
            bc_list.append(idx64)

    stdio.fclose(fin)
    stdlib.free(buffer)
    stdlib.free(line)

    dT = time() - T0
    rate = n / dT / 1000
    logging.debug(f"read {n} barcodes in {dT:.1f} seconds ({rate:.2f} k/sec)")
    return bc_list


cpdef ingest_sequences(list seq_list):
    """
    Ingest a list of byte sequences and return their uint32 encodings as a list.
    """
    cdef Py_ssize_t n = len(seq_list)
    cdef Py_ssize_t i
    cdef list result = [0] * n
    for i in range(n):
        result[i] = seq_to_uint32(seq_list[i])
    return result

@cython.cfunc
@cython.inline
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
cdef uint32_t seq_to_uint32_slice(bytes seq, Py_ssize_t start, Py_ssize_t length):
    """Compute uint32 code for seq[start:start+length] without allocating a substring."""
    cdef Py_ssize_t k, end = start + length
    cdef uint32_t result = 0
    cdef int ch
    cdef unsigned char nt
    for k in range(start, end):
        nt = seq[k]
        result <<= 2
        result |= base_map[nt]

    return result


@cython.cfunc
@cython.inline
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
cdef uint8_t find_in_list(uint32_t idx_s, uint32_t* slist, uint32_t n):
    """
    Find idx_s in the sorted suffix list slist of length n.
    Return 1 if found, 0 otherwise.
    """
    cdef uint32_t i, m
    cdef uint32_t ofs, idx_l, idx_r, idx_cur, idx_m

    if n == 0:
        return 0
    
    if n == 1:
        return slist[0] == idx_s

    idx_l = slist[0]
    idx_r = slist[n - 1]

    if idx_l == idx_s:
        return 1
    if idx_r == idx_s:
        return 1
    
    if idx_l > idx_s:
        # first entry is already larger than target
        return 0
    
    if idx_r < idx_s:
        # last entry is still smaller than target
        return 0

    if (n > 2):
        m = n // 2
        idx_m = slist[m]
        if idx_m == idx_s:
            return 1

        if (idx_m > idx_s):
            # look in left half
            return find_in_list(idx_s, &slist[1], m - 1)
        elif (idx_m < idx_s):
            # look in right half
            return find_in_list(idx_s, &slist[m + 1], n - m - 2)

    return 0


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
#@cython.nogil
def query(list bc_list, ndarray[np_uint8_t, ndim=1] hits, ndarray[np_uint32_t, ndim=1] PI, ndarray[np_uint32_t, ndim=1] SL, int l_prefix):
    """
    For each barcode in bc_list, check whether it is in the index defined by PI and SL.
    Mark hits in the hits array (1 = hit, 0 = no hit).
    bc_list: list of bytes objects (barcodes)
    hits: 1D numpy array of uint8_t, preallocated, length = len(bc_list)
    PI: 1D numpy array of uint32_t, prefix index
    SL: 1D numpy array of uint32_t, suffix list
    l_prefix: length of prefix in bases
    """
    cdef Py_ssize_t n = len(bc_list)
    cdef Py_ssize_t i
    cdef uint32_t j, idx_p, idx_s, ofs, nn
    cdef bytes bc
    # create typed memoryviews for fast C-level access
    # PI and SL can be read-only (e.g. memory-mapped files), so use const views
    cdef const uint32_t[:] PI_view = PI
    cdef const uint32_t[:] SL_view = SL
    cdef uint8_t[:] hits_view = hits

    cdef char *cbuf
    cdef Py_ssize_t blen
    # cdef const unsigned char *ubuf = <const unsigned char *> cbuf

    for i in range(n):
        bc = bc_list[i]
        # fast-path for bytes: read raw buffer and compute indices without Python indexing
        PyBytes_AsStringAndSize(bc, &cbuf, &blen)
        idx_p = seq_to_uint32_from_buf_nogil(<const unsigned char *>cbuf, 0, l_prefix)
        # print(f"idx_p={idx_p} for bc={bc}")

        ofs = PI_view[idx_p]
        # print("ofs=", ofs)
        if ofs != 0:
            idx_s = seq_to_uint32_from_buf_nogil(<const unsigned char *>cbuf, l_prefix, blen - l_prefix)
            nn = SL_view[ofs]
            # print(f"idx_s={idx_s} nn={nn}")
            if find_in_list(idx_s, &SL_view[ofs+1], nn):
                hits_view[i] = 1

            # for j in range(nn):
            #     idx_cur = SL_view[ofs + j]
            #     if idx_cur == idx_s:
            #         hits_view[i] = 1
            #         break
            #     elif idx_cur > idx_s:
            #         # indices are sorted, can break early
            #         break

    return hits





@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
#@cython.nogil
def query_idx64(list bc_list, ndarray[np_uint8_t, ndim=1] hits, ndarray[np_uint32_t, ndim=1] PI, ndarray[np_uint32_t, ndim=1] SL, int l_prefix, int l_suffix):
    """
    For each barcode in bc_list, check whether it is in the index defined by PI and SL.
    Mark hits in the hits array (1 = hit, 0 = no hit).
    bc_list: list of bytes objects (barcodes)
    hits: 1D numpy array of uint8_t, preallocated, length = len(bc_list)
    PI: 1D numpy array of uint32_t, prefix index
    SL: 1D numpy array of uint32_t, suffix list
    l_prefix: length of prefix in bases
    """
    cdef Py_ssize_t n = len(bc_list)
    cdef Py_ssize_t i
    cdef uint32_t j, idx_p, idx_s, ofs, nn
    cdef uint64_t bc
    # create typed memoryviews for fast C-level access
    # PI and SL can be read-only (e.g. memory-mapped files), so use const views
    cdef const uint32_t[:] PI_view = PI
    cdef const uint32_t[:] SL_view = SL
    cdef uint8_t[:] hits_view = hits
    cdef uint8_t rshift = 2 * l_suffix
    cdef uint64_t mask = (1 << (2 * l_suffix)) - 1

    for i in range(n):
        bc = bc_list[i]
        # fast-path for bytes: read raw buffer and compute indices without Python indexing
        idx_p = <uint32_t>(bc >> rshift)
        # print(f"idx_p={idx_p} for bc={bc}")

        ofs = PI_view[idx_p]
        # print("ofs=", ofs)

        if ofs != 0:
            idx_s = <uint32_t>(bc & mask) # lower 30 bits
            nn = SL_view[ofs]
            # print(f"idx_s={idx_s} nn={nn}")
            if find_in_list(idx_s, &SL_view[ofs+1], nn):
                hits_view[i] = 1

            # for j in range(nn):
            #     if SL_view[ofs + j + 1] == idx_s:
            #         hits_view[i] = 1
            #         break

    return hits



cdef make_shifts(uint64_t idx, uint64_t* shifts, int l=25):
    # cdef uint64_t[40] shifts
    cdef uint64_t j
    cdef int i = 0, ld_shift = 2*(l - 2), lm_shift = 2*(l-1)

    # add two-base right-shifted sequences with all possible dimers in front
    for j in range(16):
        shifts[i] = (dimers[j] << ld_shift) | (idx >> 4)
        i += 1

    # add single base right-shifted sequences
    for j in range(4):
        shifts[i] = (j << lm_shift) | (idx >> 2)
        i += 1

    # add single base left-shifted sequences
    for j in range(4):
        shifts[i] = (idx << 2) | j
        i += 1

    # now add left-shifted sequences with all possible dimers at the end
    for j in range(16):
        shifts[i] = (idx << 4) | dimers[j]
        i += 1

    # return np.array(shifts)

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
#@cython.nogil
def query_idx64_shifts(list bc_list, ndarray[np_uint8_t, ndim=1] hits, ndarray[np_uint32_t, ndim=1] PI, ndarray[np_uint32_t, ndim=1] SL, int l_prefix, int l_suffix):
    """
    For each barcode in bc_list, check whether it is in the index defined by PI and SL.
    Mark hits in the hits array (1 = hit, 0 = no hit).
    bc_list: list of bytes objects (barcodes)
    hits: 1D numpy array of uint8_t, preallocated, length = len(bc_list)
    PI: 1D numpy array of uint32_t, prefix index
    SL: 1D numpy array of uint32_t, suffix list
    l_prefix: length of prefix in bases
    """
    cdef Py_ssize_t n = len(bc_list)
    cdef Py_ssize_t i, n_shifts, k
    cdef uint32_t idx_p, idx_s, ofs, nn
    cdef uint64_t bc, bc_shift, j
    # create typed memoryviews for fast C-level access
    # PI and SL can be read-only (e.g. memory-mapped files), so use const views
    cdef const uint32_t[:] PI_view = PI
    cdef const uint32_t[:] SL_view = SL
    cdef uint8_t[:] hits_view = hits
    cdef uint64_t rshift = 2 * l_suffix
    cdef int l = l_prefix + l_suffix
    cdef uint64_t mask_l = (1 << (2 * l)) - 1
    cdef uint64_t mask = (1 << (2 * l_suffix)) - 1
    cdef uint64_t[42] shifts
    cdef int ld_shift = 2*(l - 2), lm_shift = 2*(l-1)

    for i in range(n):
        bc = bc_list[i]
        # print(f"{i} bc={bc} {uint64_to_seq(bc, l)}")

        n_shifts = 0
        # add two-base right-shifted sequences with all possible dimers in front
        for j in range(16):
            shifts[n_shifts] = (dimers[j] << ld_shift) | (bc >> 4)
            n_shifts += 1

        # add single base right-shifted sequences
        for j in range(4):
            shifts[n_shifts] = (j << lm_shift) | (bc >> 2)
            n_shifts += 1

        # add single base left-shifted sequences
        for j in range(4):
            shifts[n_shifts] = ((bc << 2) & mask_l) | j
            n_shifts += 1

        # now add left-shifted sequences with all possible dimers at the end
        for j in range(16):
            shifts[n_shifts] = ((bc << 4) & mask_l) | dimers[j]
            n_shifts += 1

        # print(f"made {n_shifts}")
        # make_shifts(bc, &shifts, l_prefix)
        # fast-path for bytes: read raw buffer and compute indices without Python indexing
        
        n_hits = 0
        for k in range(n_shifts):
            bc_shift = shifts[k]
            # print(f"bc_shift {k}={bc_shift} {uint64_to_seq(bc_shift, l)}")
            idx_p = <uint32_t>(bc_shift >> rshift)
            # print(f"idx_p={idx_p} -> seq={uint64_to_seq(idx_p, l_prefix)}")

            ofs = PI_view[idx_p] # Could this have side-effects if we modify ofs?
            # print(f"{i}: query {k}/{n_shifts} -> ofs={ofs}")
            # print("ofs=", ofs)
            if ofs != 0:
                idx_s = <uint32_t>(bc_shift & mask) # lower 30 bits
                nn = SL_view[ofs]
                # print("nn=", nn)
                for j in range(nn):
                    if SL_view[ofs + j + 1] == idx_s:
                        n_hits += 1
                        break

        hits_view[i] = n_hits
        # print(f"{i} -> n_hits={n_hits}")

    return hits

@cython.cfunc
@cython.inline
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
cdef make_variants_off_by_one(uint32_t idx_s, uint32_t* variants, int l=15):
    """
    Generate all single-base off-by-one variants of a suffix index.
    Return as an array of 75 uint32_t values.
    """
    
    cdef int pos, j, shift
    cdef uint32_t mut_code, idx_v

    cdef int i = 0
    for pos in range(l):  # up to 15 bases in suffix
        shift = 2 * pos
        idx_v = idx_s
        for j in range(3):
            mut_code = ((idx_v >> shift) + 1) & 0b11
            idx_v = idx_s & ~(0b11 << shift) | (mut_code << shift)
            variants[i] = idx_v
            i += 1

@cython.cfunc
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
@cython.inline
cdef make_insertions(uint64_t idx, uint64_t* variants, uint64_t l=25):
    """
    Generate all single-base insertions/clipped to l again.
    Write into an array of l*4 uint32_t values.
    """
    
    cdef uint64_t idx_l, idx0, idx_var, k, mask_l, mask_r, pos
    cdef uint64_t full_mask = 1
    full_mask = full_mask << 2*l
    full_mask -= 1

    # print(f"{uint64_to_seq(idx, l)} one base insertions")
    # print(f"l={l} 1 << 2*l = {1 << (2*l):b}")
    # print(f"full_mask {full_mask:b}")

    cdef int i = 0
    idx_l = (idx << 2) & full_mask # on base left-shifted and clipped
    idx_r = idx

    mask_l = full_mask ^ 3
    mask_r = 0
    # print(f"{uint64_to_seq(idx_l, l)} left-shifted index")
    for pos in range(l):
        # insertion at base pos
        idx0 = (idx_l & mask_l) | (idx & mask_r)
        
        # print(f"{mask_l:050b} mask_l")
        # print(f"{mask_r:050b} mask_r")
        # print(f"{uint64_to_seq(idx0, l)} IDX0 pos={pos} i={i}")

        for k in range(4):
            idx_var = idx0 | (k << pos*2)
            variants[i] = idx_var
            i += 1
        
            # print(f"{uint64_to_seq(idx_var, l)} k={k} i={i}")

        mask_l = (mask_l << 2) & full_mask
        mask_r = (mask_r << 2) | 3

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
#@cython.nogil
def query_idx64_off_by_one(list bc_list, ndarray[np_uint8_t, ndim=1] hits, ndarray[np_uint32_t, ndim=1] PI, ndarray[np_uint32_t, ndim=1] SL, int l_prefix, int l_suffix):
    """
    For each barcode in bc_list, check whether it is in the index defined by PI and SL.
    Mark hits in the hits array (1 = hit, 0 = no hit).
    bc_list: list of bytes objects (barcodes)
    hits: 1D numpy array of uint8_t, preallocated, length = len(bc_list)
    PI: 1D numpy array of uint32_t, prefix index
    SL: 1D numpy array of uint32_t, suffix list
    l_prefix: length of prefix in bases
    """
    cdef Py_ssize_t n = len(bc_list)
    cdef Py_ssize_t i, j, k, n_p_variants = l_prefix * 3, n_s_variants = l_suffix * 3
    cdef uint32_t idx_p, idx_s, idx_sv, ofs, nn
    cdef uint64_t bc
    # create typed memoryviews for fast C-level access
    # PI and SL can be read-only (e.g. memory-mapped files), so use const views
    cdef const uint32_t[:] PI_view = PI
    cdef const uint32_t[:] SL_view = SL
    cdef uint8_t[:] hits_view = hits
    cdef uint8_t rshift = 2 * l_suffix
    cdef uint64_t mask = (1 << (2 * l_suffix)) - 1
    cdef uint32_t[16*3] idx_s_variants # 15 * 3 = 45 variants for single-base off-by-one
    cdef uint32_t[16*3] idx_p_variants

    for i in range(n):
        bc = bc_list[i]

        # fast-path for bytes: read raw buffer and compute indices without Python indexing
        idx_p = <uint32_t>(bc >> rshift)

        ofs = PI_view[idx_p]
        n_hits = 0
        if ofs != 0:
            idx_s = <uint32_t>(bc & mask) # lower 2 x l_suffix bits
            make_variants_off_by_one(idx_s, idx_s_variants, l_suffix)
            nn = SL_view[ofs]
            ofs += 1
            for idx_sv in idx_s_variants[:n_s_variants]:
                if find_in_list(idx_sv, &SL_view[ofs+1], nn):
                    n_hits += 1
        
        # now also check prefix variants
        make_variants_off_by_one(idx_p, idx_p_variants, l_prefix)
        for idx_p in idx_p_variants[:n_p_variants]:
            ofs = PI_view[idx_p]
            if ofs != 0:
                nn = SL_view[ofs]
                if find_in_list(idx_s, &SL_view[ofs+1], nn):
                    n_hits += 1

        hits_view[i] = n_hits   

    return hits



@cython.wraparound(False)
@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
@cython.exceptval(check=False)
# @cython.nogil
def query_idx64_indel(list bc_list, ndarray[np_uint8_t, ndim=1] hits, ndarray[np_uint32_t, ndim=1] PI, ndarray[np_uint32_t, ndim=1] SL, int l_prefix, int l_suffix):
    """
    For each barcode in bc_list, check whether it is in the index defined by PI and SL.
    Mark hits in the hits array (1 = hit, 0 = no hit).
    bc_list: list of bytes objects (barcodes)
    hits: 1D numpy array of uint8_t, preallocated, length = len(bc_list)
    PI: 1D numpy array of uint32_t, prefix index
    SL: 1D numpy array of uint32_t, suffix list
    l_prefix: length of prefix in bases
    """
    cdef Py_ssize_t n = len(bc_list)
    cdef Py_ssize_t i
    cdef uint32_t j, idx_p, idx_s, ofs, nn
    cdef uint64_t bc, bc_var, l = l_prefix + l_suffix
    # create typed memoryviews for fast C-level access
    # PI and SL can be read-only (e.g. memory-mapped files), so use const views
    cdef const uint32_t[:] PI_view = PI
    cdef const uint32_t[:] SL_view = SL
    cdef uint8_t[:] hits_view = hits
    cdef uint8_t rshift = 2 * l_suffix
    cdef uint64_t mask = (1 << (2 * l_suffix)) - 1, n_indels = 4*l # for now just the insertions

    cdef uint64_t[32*4] idx_indels

    for i in range(n):
        bc = bc_list[i]
        make_insertions(bc, idx_indels, l)

        for bc_var in idx_indels[:n_indels]:
            # fast-path for bytes: read raw buffer and compute indices without Python indexing
            idx_p = <uint32_t>(bc_var >> rshift)
            # print(f"idx_p={idx_p} for bc={bc_var}")

            ofs = PI_view[idx_p]
            # print("ofs=", ofs)

            if ofs != 0:
                idx_s = <uint32_t>(bc_var & mask) # lower 30 bits
                nn = SL_view[ofs]
                # print(f"idx_s={idx_s} nn={nn}")
                if find_in_list(idx_s, &SL_view[ofs+1], nn):
                    hits_view[i] += 1

                # for j in range(nn):
                #     if SL_view[ofs + j + 1] == idx_s:
                #         hits_view[i] = 1
                #         break

    return hits