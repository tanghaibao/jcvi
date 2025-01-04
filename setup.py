#!/usr/bin/env python

"""Package setup for Cython extensions only"""

from Cython.Build import build_ext
from setuptools import setup, Extension
import numpy as np

ext_modules = [
    Extension(
        "jcvi.assembly.chic",
        ["src/jcvi/assembly/chic.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3"],
    ),
    Extension(
        "jcvi.formats.cblast",
        ["src/jcvi/formats/cblast.pyx"],
        extra_compile_args=["-O3"],
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
