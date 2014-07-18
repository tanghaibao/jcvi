#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
      name='jcvi',
      version='0.4.7',
      author='Haibao Tang',
      author_email='tanghaibao@gmail.com',
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False,
      license='LICENSE',
      url='http://pypi.python.org/pypi/jcvi/',
      description='Python utility libraries on genome assembly, '\
                  'annotation and comparative genomics',
      long_description=open("README.rst").read(),
      install_requires=['biopython', 'numpy', 'matplotlib']
     )
