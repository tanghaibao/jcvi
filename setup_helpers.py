#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import subprocess
import pkg_resources


def check_version(majorv=2, minorv=7):
    """ Make sure the package runs on the supported Python version
    """
    if sys.version_info.major == majorv and sys.version_info.minor != minorv:
        sys.stderr.write("ERROR: %s is only for >= Python %d.%d but you are running %d.%d\n" %\
                    (name, majorv, minorv, sys.version_info.major, sys.version_info.minor))
        sys.exit(1)


def get_init(filename="__init__.py"):
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


def get_long_description(filename='README.md'):
    """ I really prefer Markdown to reStructuredText. PyPi does not.
    """
    try:
       import pypandoc
       description = pypandoc.convert('README.md', 'rst')
    except (IOError, ImportError):
       description = open("README.md").read()
    return description
