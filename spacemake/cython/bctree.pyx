#cython: boundscheck=False, wraparound=False, initializedcheck=False, overflowcheck=False, cdivision=True, language_level=3
##cython: boundscheck=True, wraparound=True, initializedcheck=True, overflowcheck=True, cdivision=False, language_level=3
#!python

__license__ = "MIT"
__authors__ = ["Marvin Jens"]
__email__ = "marvin.jens@charite.de"

from types cimport *
from types import *

from cython.parallel import parallel, prange
import numpy as np
cimport numpy as np
cimport cython
cimport openmp

from libc.math cimport exp, log
from libc.stdlib cimport abort, malloc, free

# maps ASCII values of A,C,G,T (U) to correct bits
# N is silently converted to A
cdef UINT8_t letter_to_bits[256]
for i in range(256):
    letter_to_bits[i] = 255

letter_to_bits[ord('n')] = 0
letter_to_bits[ord('a')] = 0 
letter_to_bits[ord('c')] = 1 
letter_to_bits[ord('g')] = 2
letter_to_bits[ord('t')] = 3
letter_to_bits[ord('u')] = 3

letter_to_bits[ord('N')] = 0
letter_to_bits[ord('A')] = 0 
letter_to_bits[ord('C')] = 1 
letter_to_bits[ord('G')] = 2
letter_to_bits[ord('T')] = 3
letter_to_bits[ord('U')] = 3


cdef UINT8_t bits_to_letters[4]
bits_to_letters[:] = [ord('A'), ord('C'), ord('G'), ord('T')]
    
cpdef inline seq_to_index(str py_str):
    cdef UINT32_t L = len(py_str)
    cdef UINT8_t x = 0
    cdef UINT8_t nt = 0
    #cdef bytes py_byte_str = 
    cdef const unsigned char[:] seq = py_str.encode('ascii')
    
    cdef UINT64_t idx = 0
    for x in range(L):
        idx = (idx << 2) | letter_to_bits[seq[x]]

    return idx

def index_to_seq(UINT64_t index, UINT8_t l):
    nucs = ['A','C','G','T']
    seq = []
    for i in range(l):
        j = index >> ((l-i-1) * 2)
        seq.append(nucs[j & 3])

    return "".join(seq)

cpdef np.ndarray[UINT64_t] hull(UINT64_t index, UINT8_t l):
    """
    Generate all 3 x l edit-distance 1 variants of the 2-bit encoded l-mer index
    """
    cdef np.ndarray[UINT64_t] variants = np.zeros(3*l + 1, dtype=np.uint64)
    cdef UINT64_t i, j, var, mask, mut, base
   
    # first entry is the unedited reference
    variants[0] = index

    for i in range(l):
        # all bits 1 except the two bits of the base we are editing
        mask = ~(<unsigned long> 3 << (i << 1))
        # NOTE: i << 1 is the same as i * 2
        for j in range(3):
            # select the bits NOT to edit
            var = index & mask
            # select the bits of the base we DO want to edit (and bring to lower bits)
            base = (index >> (i << 1)) & <unsigned long> 3
            # addition plus wrap around at 3 gives all mutations
            mut = (base + 1 + j) % <unsigned long> 4 
            # add in the edited base bits at the correct bit position again
            var |= mut << (i * 2)
            # done. store the edited l-mer
            variants[i * <unsigned long> 3 + j + 1] = var
    
    return variants

cdef packed struct bc_info_t:
    # Data structure for the barcode tree
    UINT64_t idx # can store barcodes up to 32nt length     # 8
    UINT32_t left # can store up to 4 billion barcodes    # 4
    UINT32_t right                                        # 4
    UINT16_t tile # up to 65k tiles                     # 2
    UINT16_t x # x/y positions between 0..65k           # 2
    UINT16_t y                                          # 2
    # total size = 22 bytes

ctypedef bc_info_t bc_info

cdef UINT32_t LEAF = 0xFFFFFFFE

cpdef UINT32_t read_txt_to_buf(src, bc_info_t [:] buf, UINT16_t tile=0):
    #cdef bc_info_t [::1] buf = _buf # fast memoryview into the buffer
    cdef str line
    cdef list cols
    #cdef bc_info_t[] records = malloc(sizeof(bc_info_t) * N_max)
    cdef bc_info_t* bc
    cdef UINT64_t i = 0
    for line in src:
        if line.startswith('cell_bc'):
            continue
        
        cols = line.split('\t')
        bc = &buf[i]
        bc.idx = seq_to_index(cols[0])
        bc.tile = tile
        bc.x = int(cols[1])
        bc.y = int(cols[2])
        bc.left = LEAF
        bc.right = LEAF

        # print(i, buf[i])

        i += 1
    
    return i

cpdef UINT32_t build_tree(bc_info_t [:] buf, UINT32_t low, UINT32_t high):
    cdef bc_info_t* node
    cdef UINT32_t root = 0

    if high - low < 1:
        return LEAF

    cdef UINT32_t mid = ((high - low) >> 1) + low
    node = &buf[mid]
    node.left = build_tree(buf, low, mid)
    node.right = build_tree(buf, mid + 1, high)
   
    return mid

cpdef inline UINT8_t check_contains(const bc_info_t [:] buf, UINT32_t root, UINT64_t idx):
    cdef bc_info_t* node
    cdef bc_info_t* leaf
    if root == LEAF:
        return 0
    else:
        node = &buf[root]
        if node.idx == idx:
            return True
        elif node.idx > idx:
            return check_contains(buf, node.left, idx)
        elif node.idx < idx:
            return check_contains(buf, node.right, idx)

cpdef np.ndarray[UINT8_t] check_contains_batch(const bc_info_t [:] buf, UINT32_t root, const UINT64_t [:] indices):
    cdef UINT64_t N = len(indices)

    # output buffer
    cdef np.ndarray[UINT8_t] res_ = np.empty(N, dtype=np.uint8)
    cdef UINT8_t [::1] res = res_ # fast memory-view object
    cdef UINT64_t i
    for i in range(N):
        res[i] = check_contains(buf, root, indices[i])

    return res_


