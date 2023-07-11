__license__ = "MIT"
__version__ = "0.9.8"
__authors__ = ["Marvin Jens"]
__email__ = "mjens@mit.edu"


# from cython.parallel import parallel, prange
import numpy as np
cimport numpy as np
# cimport cython
# cimport openmp

# from libc.math cimport exp, log
# from libc.stdlib cimport abort, malloc, free

ctypedef extern np.uint8_t UINT8_t
ctypedef extern np.uint16_t UINT16_t
ctypedef extern np.uint32_t UINT32_t
ctypedef extern np.uint64_t UINT64_t
ctypedef extern np.int32_t INT32_t
ctypedef extern np.float32_t FLOAT32_t
ctypedef extern np.float64_t FLOAT64_t
                                
# cdef extern UINT64_t rand_state[2]
# cdef extern UINT64_t RAND_MAX = 2**64 - 1
# cdef extern FLOAT32_t FRAND_MAX = RAND_MAX

# cdef extern UINT8_t letter_to_bits[256]