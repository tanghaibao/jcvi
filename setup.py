#!/usr/bin/env python

"""Package setup"""

from Cython.Build import build_ext
from setuptools import setup, Extension

import numpy as np

ext_modules = [
    Extension(
        "jcvi.assembly.chic",
        ["jcvi/assembly/chic.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3"],
    ),
    Extension(
        "jcvi.formats.cblast", ["jcvi/formats/cblast.pyx"], extra_compile_args=["-O3"]
    ),
]

if __name__ == "__main__":
    setup(
        cmdclass={"build_ext": build_ext},
        ext_modules=ext_modules,
    )
