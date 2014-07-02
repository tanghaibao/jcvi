#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='jcvi',
      version='1.0',
      description='Python utility libraries on genome assembly, '\
                  'annotation and comparative genomics',
      author='Haibao Tang',
      author_email='tanghaibao@gmail.com',
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False,
     )
