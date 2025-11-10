# cython: boundscheck=False, wraparound=False, cdivision=True, nonecheck=False, initializedcheck=False
## cython: language_level=3
import cython
cimport cython
from libc.stdint cimport uint32_t, uint64_t
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


import numpy as np
from libc.stdint cimport uint8_t
from numpy cimport ndarray, uint32_t as np_uint32_t, uint8_t as np_uint8_t

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

        ofs = PI_view[idx_p]
        if ofs != 0:
            idx_s = seq_to_uint32_from_buf_nogil(<const unsigned char *>cbuf, l_prefix, blen - l_prefix)
            nn = SL_view[ofs]
            ofs += 1
            for j in range(nn):
                if SL_view[ofs + j] == idx_s:
                    hits_view[i] = 1
                    break

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

        ofs = PI_view[idx_p]
        if ofs != 0:
            idx_s = <uint32_t>(bc & mask) # lower 30 bits
            nn = SL_view[ofs]
            ofs += 1
            for j in range(nn):
                if SL_view[ofs + j] == idx_s:
                    hits_view[i] = 1
                    break

    return hits
