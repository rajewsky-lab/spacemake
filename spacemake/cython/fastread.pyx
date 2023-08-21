#cython: boundscheck=False, wraparound=False, initializedcheck=False, overflowcheck=False, cdivision=True, language_level=3
##cython: boundscheck=True, wraparound=True, initializedcheck=True, overflowcheck=True, cdivision=False, language_level=3
#!python

from types cimport *
from types import *

from cython.parallel import parallel, prange
import numpy as np
cimport numpy as np
cimport cython
cimport openmp

from libc.math cimport exp, log
from libc.stdlib cimport abort, malloc, free

def read_txt_chunked(src, UINT64_t chunk_size=100000):
    cdef str line
    cdef list buf = []
    cdef UINT64_t i = 0
    for line in src:
        buf.append(line)
        i += 1
        if i >= chunk_size:
            yield buf
            buf = []
            i = 0

    if i > 0:
        yield buf
        buf = []
        i = 0
