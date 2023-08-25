#cython: boundscheck=False, wraparound=False, initializedcheck=False, overflowcheck=False, cdivision=True
###cython: boundscheck=True, wraparound=True, initializedcheck=True, overflowcheck=True, cdivision=False
#!python

#from types cimport *
from types import *
from libc cimport stdlib, stdio
from typing import TextIO

def distribute(str fin_name, list fifo_names, int chunk_size=10000):
    # int fd_input = stdlib.open(finput, mode=stdlib.O_RDONLY)
    # wh
    cdef int i = 0
    cdef int j = 0
    cdef size_t n_outs = len(fifo_names)
    cdef size_t n = 0
    cdef str line
    cdef char* buffer
    cdef size_t L = 2**20 # 512 kb
    cdef size_t L_stdin = 2**20 # 1Mb
    cdef ssize_t n_read
    buffer = <char*>stdlib.malloc(L)
    cdef char* fifo_buffers[32]

    assert n <= 32
    cdef stdio.FILE *fifos[32]

    cdef char* stdin_buf = <char*>stdlib.malloc(L_stdin)
    cdef stdio.FILE *fin = stdio.fopen(fin_name.encode('utf-8'), 'r')
    stdio.setvbuf(fin, stdin_buf, stdio._IOFBF, L_stdin)

    for i in range(n_outs):
        fifo_buffers[i] = <char*>stdlib.malloc(L)
        fifos[i] = stdio.fopen(fifo_names[i].encode('utf-8'), 'w')
        stdio.setvbuf(fifos[i], fifo_buffers[i], stdio._IOFBF, L)

    while(True):
        n_read = stdio.getline(&buffer, &L, fin) 
        if n_read <= 0:
            break

        j = n // chunk_size
        stdio.fwrite(buffer, n_read, 1, fifos[j % n_outs])
        n += 1

    for i in range(n_outs):
        stdio.fclose(fifos[i])
        stdlib.free(fifo_buffers[i])

    stdlib.free(buffer)
    print(n)


