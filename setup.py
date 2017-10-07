#!/usr/bin/env python

import sys
import subprocess
import pkg_resources

from setuptools import setup, find_packages, Extension
from setup_helpers import check_version, get_init, missing_requirements, \
            install_requirements, get_long_description


name = "jcvi"
classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]


# Temporarily install dependencies required by setup.py before trying to import them.
sys.path[0:0] = ['setup-requires']
pkg_resources.working_set.add_entry('setup-requires')

requires = ["cython", "numpy"]
install_requirements(missing_requirements(requires))


# Now these are available
import numpy as np
from Cython.Distutils import build_ext

# Start the show
check_version()
author, email, license, version = get_init()
ext_modules = [
    Extension("jcvi.assembly.chic", ["assembly/chic.pyx"],
                     include_dirs=[np.get_include()],
                     extra_compile_args=["-O3"]),
    Extension("jcvi.formats.cblast", ["formats/cblast.pyx"],
                     extra_compile_args=["-O3"])
]


setup(
      name=name,
      author=author[0],
      author_email=email,
      version=version,
      cmdclass={'build_ext': build_ext},
      package_dir={name: '.'},
      packages=[name] + ['.'.join((name, x)) for x in find_packages()],
      include_package_data=True,
      package_data={"jcvi.utils.data": ["*.*"]},
      ext_modules=ext_modules,
      classifiers=classifiers,
      zip_safe=False,
      license='BSD',
      url='http://github.com/tanghaibao/jcvi',
      description='Python utility libraries on genome assembly, '\
                  'annotation and comparative genomics',
      long_description=get_long_description(),
      install_requires=['biopython', 'deap',
                        'matplotlib', 'networkx', 'numpy'],
 )
