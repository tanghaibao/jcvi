#!/usr/bin/env python

import sys
import subprocess
import pkg_resources

from setuptools import setup, find_packages, Extension


name = "jcvi"
majorv, minorv = 2, 7
classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]


def check_version():
    """ Make sure the package runs on the supported Python version
    """
    if sys.version_info.major == majorv and sys.version_info.minor != minorv:
        sys.stderr.write("ERROR: %s is only for >= Python %d.%d but you are running %d.%d\n" %\
                    (name, majorv, minorv, sys.version_info.major, sys.version_info.minor))
        sys.exit(1)


def import_init(filename="__init__.py"):
    """ Get various info from the package without importing them
    """
    import ast

    with open(filename) as init_file:
        module = ast.parse(init_file.read())

    itr = lambda x: (ast.literal_eval(node.value) for node in ast.walk(module) \
        if isinstance(node, ast.Assign) and node.targets[0].id == x)

    try:
        return next(itr("__author__")), \
               next(itr("__email__")), \
               next(itr("__license__")), \
               next(itr("__version__"))
    except StopIteration:
        raise ValueError("One of author, email, license, or version"
                    " cannot be found in {}".format(filename))


def missing_requirements(specifiers):
    """ Find what's missing
    """
    for specifier in specifiers:
        try:
            pkg_resources.require(specifier)
        except pkg_resources.DistributionNotFound:
            yield specifier


def install_requirements(specifiers):
    """ Install the listed requirements
    """
    to_install = list(specifiers)
    if to_install:
        cmd = [sys.executable, "-m", "pip", "install",
            "-t", "setup-requires"] + to_install
        subprocess.call(cmd)


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
author, email, license, version = import_init()
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
      long_description=open("README.md").read(),
      install_requires=['biopython', 'deap',
                        'matplotlib', 'networkx', 'numpy'],
 )
