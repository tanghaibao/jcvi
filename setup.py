#!/usr/bin/env python

from __future__ import absolute_import

import os.path as op
import versioneer

from setuptools import setup, find_packages, Extension
from setup_helper import SetupHelper

name = "jcvi"
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

# Use the helper
h = SetupHelper(initfile="jcvi/__init__.py", readmefile="README.md")
h.check_version(name, majorv=2, minorv=7)
cmdclass = versioneer.get_cmdclass()
include_dirs = []
setup_dir = op.abspath(op.dirname(__file__))
requirements = [x.strip() for x in open(op.join(setup_dir, "requirements.txt"))]
h.install_requirements(requires=["cython", "numpy"])

# Build the ext
try:
    import numpy as np
    from Cython.Distutils import build_ext

    cmdclass.update({"build_ext": build_ext})
    include_dirs.append(np.get_include())
except ImportError:
    print("Cython not installed. Skip compiling Cython extensions.")

ext_modules = [
    Extension(
        "jcvi.assembly.chic",
        ["jcvi/assembly/chic.pyx"],
        include_dirs=include_dirs,
        extra_compile_args=["-O3"],
    ),
    Extension(
        "jcvi.formats.cblast", ["jcvi/formats/cblast.pyx"], extra_compile_args=["-O3"]
    ),
]

packages = [name] + [
    ".".join((name, x)) for x in find_packages("jcvi", exclude=["test*.py"])
]

setup(
    name=name,
    author=h.author,
    author_email=h.email,
    version=versioneer.get_version(),
    license=h.license,
    long_description=h.long_description,
    long_description_content_type="text/markdown",
    cmdclass=cmdclass,
    packages=packages,
    include_package_data=True,
    package_data={"jcvi.utils.data": ["*.*"]},
    ext_modules=ext_modules,
    classifiers=classifiers,
    zip_safe=False,
    url="http://github.com/tanghaibao/jcvi",
    description="Python utility libraries on genome assembly, "
    "annotation and comparative genomics",
    setup_requires=["setuptools>=18.0", "cython"],
    install_requires=requirements,
)
