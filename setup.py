#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
from setup_helper import SetupHelper


name = "jcvi"
classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]

# Use the helper
h = SetupHelper(initfile="__init__.py", readmefile="README.md")
h.check_version(name, majorv=2, minorv=7)
h.install_requirements(requires=["cython", "numpy"])

# Now these are available
import numpy as np
from Cython.Distutils import build_ext

# Start the show
ext_modules = [
    Extension("jcvi.assembly.chic", ["assembly/chic.pyx"],
                     include_dirs=[np.get_include()],
                     extra_compile_args=["-O3"]),
    Extension("jcvi.formats.cblast", ["formats/cblast.pyx"],
                     extra_compile_args=["-O3"])
]

setup(
      name=name,
      author=h.author,
      author_email=h.email,
      version=h.version,
      license=h.license,
      long_description=h.long_description,
      cmdclass={'build_ext': build_ext},
      package_dir={name: '.'},
      packages=[name] + ['.'.join((name, x)) for x in find_packages()],
      include_package_data=True,
      package_data={"jcvi.utils.data": ["*.*"]},
      ext_modules=ext_modules,
      classifiers=classifiers,
      zip_safe=False,
      url='http://github.com/tanghaibao/jcvi',
      description='Python utility libraries on genome assembly, '\
                  'annotation and comparative genomics',
      install_requires=['biopython', 'deap',
                        'matplotlib', 'networkx', 'numpy'],
 )
