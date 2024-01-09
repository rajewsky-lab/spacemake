from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Compiler import Options
import numpy

# These are optional
Options.docstrings = True
Options.annotate = False

# Modules to be compiled and include_dirs when necessary
extensions = [
    # Extension(
    #     "pyctmctree.inpyranoid_c",
    #     ["src/pyctmctree/inpyranoid_c.pyx"],
    # ),
    Extension(
        "spacemake.cython.fast_loop",
        ["spacemake/cython/fast_loop.pyx"], include_dirs=[numpy.get_include()],
    ),
    Extension(
        "spacemake.cython.fastread",
        ["spacemake/cython/fastread.pyx"], include_dirs=[numpy.get_include()],
    ),

]

if __name__ == "__main__":
    from setuptools import setup

    setup(
        name='spacemake',  # Required

        # A list of compiler Directives is available at
        # https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#compiler-directives

        # external to be compiled
        ext_modules = cythonize(extensions, compiler_directives={"language_level": 3, "profile": False}),
    )
