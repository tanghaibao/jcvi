#!/usr/bin/env python

from setuptools import setup, find_packages
from glob import glob


name = "jcvi"
classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]

exec(open("__init__.py").read())
setup(
      name=name,
      version=__version__,
      author=__author__[0],
      author_email=__email__,
      package_dir={name: '.'},
      packages=[name] + ['.'.join((name, x)) for x in find_packages()],
      include_package_data=True,
      package_data={"jcvi.utils.data": ["*.*"]},
      classifiers=classifiers,
      zip_safe=False,
      license='BSD',
      url='http://github.com/tanghaibao/jcvi',
      description='Python utility libraries on genome assembly, '\
                  'annotation and comparative genomics',
      long_description=open("README.rst").read(),
      install_requires=['biopython', 'numpy', 'matplotlib',
                        'deap', 'networkx']
     )
