#cython: boundscheck=False, wraparound=False, initializedcheck=False, overflowcheck=False, cdivision=True
###cython: boundscheck=True, wraparound=True, initializedcheck=True, overflowcheck=True, cdivision=False
#!python

#from types cimport *
from types import *
from libc cimport stdlib, stdio
from libc.string cimport memcpy
from typing import TextIO

def distribute(str fin_name, list fifo_names, int chunk_size=10000, size_t in_buf_size=2**20, size_t out_buf_size=2**19, header_detect_func=None, header_fifo=None):
    cdef size_t i = 0
    cdef size_t j = 0
    cdef size_t n_outs = len(fifo_names)
    cdef size_t n = 0
    cdef str line
    cdef char* buffer
    cdef ssize_t n_read = 0
    # support input distribution to up to 128 fifos in parallel
    cdef char* fifo_buffers[128]
    buffer = <char*>stdlib.malloc(in_buf_size)

    assert n <= 128
    cdef stdio.FILE *fifos[128]

    cdef char* stdin_buf = <char*>stdlib.malloc(in_buf_size)
    cdef stdio.FILE *fin = stdio.fopen(fin_name.encode('utf-8'), 'r')
    cdef stdio.FILE *fheader
    stdio.setvbuf(fin, stdin_buf, stdio._IOFBF, in_buf_size)


    if header_detect_func and header_fifo:
        # if the input stream contains some special header lines, allow to detect these and write them to a separate header-fifo
        fheader = stdio.fopen(header_fifo.encode('utf-8'), 'w')
        while(True):
            n_read = stdio.getline(&buffer, &in_buf_size, fin) 
            if n_read <= 0:
                break
            
            if header_detect_func(str(buffer)):
                stdio.fwrite(buffer, n_read, 1, fheader)
            else:
                # as soon as we see a non-header line, we break and enter the main loop
                break

        # close the header fifo
        stdio.fclose(fheader)

    # set up output buffers for the main fifos
    for i in range(n_outs):
        fifo_buffers[i] = <char*>stdlib.malloc(out_buf_size)
        fifos[i] = stdio.fopen(fifo_names[i].encode('utf-8'), 'w')
        stdio.setvbuf(fifos[i], fifo_buffers[i], stdio._IOFBF, out_buf_size)

    # main distribute loop
    while(True):
        if n_read == 0:
            # we may have already read a line in the header-detection above!
            # read a new-line only if n_read has explicitly been zeroed again
            n_read = stdio.getline(&buffer, &in_buf_size, fin) 

        if n_read <= 0:
            break

        j = n // chunk_size
        stdio.fwrite(buffer, n_read, 1, fifos[j % n_outs])
        n += 1
        n_read = 0

    # close down all main fifos and free the output buffers
    for i in range(n_outs):
        stdio.fclose(fifos[i])
        stdlib.free(fifo_buffers[i])

    # close the input and free the input buffer, as well as the line buffer
    stdio.fclose(fin)
    stdlib.free(stdin_buf)
    stdlib.free(buffer)

    return n # number of lines distributed (excluding header)


def collect(list fifo_names, str fout_name, int chunk_size=10000, size_t in_buf_size=2**19, size_t out_buf_size=2**20, header_fifo=None):
    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    cdef size_t n_ins = len(fifo_names)
    cdef size_t n = 0
    cdef str line
    cdef char* buffer
    cdef ssize_t n_read
    
    # support collecting and demuxing from up to 128 fifos
    cdef char* fifo_buffers[128]

    import sys
    # sys.stderr.write(f"collecting from {fifo_names} into {fout_name}")
    assert n <= 128
    cdef stdio.FILE *fifos[128]

    
    # line buffer
    buffer = <char*>stdlib.malloc(out_buf_size)
    # prepare output file and buffer
    cdef char* out_buf = <char*>stdlib.malloc(out_buf_size)
    cdef stdio.FILE *fout = stdio.fopen(fout_name.encode('utf-8'), 'w')
    stdio.setvbuf(fout, out_buf, stdio._IOFBF, out_buf_size)

    cdef stdio.FILE *fheader
    if header_fifo:
        # we have some special, un-multiplexed header data. Read this stuff first and write it to the output
        fheader = stdio.fopen(header_fifo.encode('utf-8'), 'r')
        while(True):
            n_read = stdio.getline(&buffer, &out_buf_size, fheader) 
            if n_read <= 0:
                break
            
            stdio.fwrite(buffer, n_read, 1, fout)
        
        # all header data has been read! close this fifo's reading end
        stdio.fclose(fheader)

    # now prepare input from all the demuxed fifos
    for i in range(n_ins):
        fifo_buffers[i] = <char*>stdlib.malloc(in_buf_size)
        fifos[i] = stdio.fopen(fifo_names[i].encode('utf-8'), 'r')
        stdio.setvbuf(fifos[i], fifo_buffers[i], stdio._IOFBF, in_buf_size)

    # main mutiplexing loop
    while(True):
        j = n // chunk_size
        k = j % n_ins # index of fifo
        # print(f"collection wants to read from {fifo_names[k]}. blocking")
        n_read = stdio.getline(&buffer, &out_buf_size, fifos[k]) 
        if n_read <= 0:
            break
        
        # print(f"fifo {k} ({fifo_names[k]}). {n_read} bytes. writing to {fout_name}")
        stdio.fwrite(buffer, n_read, 1, fout)
        n += 1

    for i in range(n_ins):
        stdio.fclose(fifos[i])
        stdlib.free(fifo_buffers[i])

    # free the line buffer, close the output, and free the output buffer
    stdlib.free(buffer)
    stdio.fclose(fout)
    stdlib.free(out_buf)

    return n # total number of multiplexed lines (excluding header)
