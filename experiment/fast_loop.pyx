#cython: boundscheck=False, wraparound=False, initializedcheck=False, overflowcheck=False, cdivision=True
###cython: boundscheck=True, wraparound=True, initializedcheck=True, overflowcheck=True, cdivision=False
#!python

#from types cimport *
from types import *
from libc cimport stdlib, stdio

def distribute(finput, fifos, int rec_size=4, int chunk_size=10000):
    # int fd_input = stdlib.open(finput, mode=stdlib.O_RDONLY)
    # wh
    cdef int i = 0
    cdef int j = 0
    cdef int n = len(fifos)
    cdef str line

    for line in finput:
        j = i % n
        fifos[j].write(line)
        # stdio.write(fifos[j].fileno, line)
        i += 1



