from setuptools import setup, Extension
import numpy as np
from Cython.Build import cythonize

extensions = [
            Extension(
                        "seqidx", ["seqidx.pyx"], include_dirs=[np.get_include(),], extra_compile_args = ["-O3", "-march=native"],
                      )
            ]

setup(
            name="spacemake",
            ext_modules=cythonize(extensions),
            zip_safe=False,
)
#import setuptools
#import Cython.Build as cb
#setuptools.setup(name='My Cython Project',     ext_modules=cb.cythonize('seqidx.pyx'))
