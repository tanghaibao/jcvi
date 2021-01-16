#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""Helper functions for package setup.
"""

import ast
import sys
import subprocess
import pkg_resources


class SetupHelper(object):
    """
    SetupHelper streamlines the population of setup() function calls with
    contents from __init__.py and README.md.
    """

    def __init__(self, initfile="__init__.py", readmefile="README.md"):
        self.author, self.email, self.license = self.get_init(initfile)
        self.author = ", ".join(self.author)
        self.long_description = self.get_long_description(readmefile)

    def check_version(self, name, majorv=2, minorv=7):
        """Make sure the package runs on the supported Python version"""
        if sys.version_info < (majorv, minorv):
            sys.stderr.write(
                "ERROR: %s is only for >= Python %d.%d but you are running %d.%d\n"
                % (name, majorv, minorv, sys.version_info.major, sys.version_info.minor)
            )
            sys.exit(1)

    def get_init(self, filename="__init__.py"):
        """Get various info from the package without importing them"""
        with open(filename) as init_file:
            module = ast.parse(init_file.read())

        itr = lambda x: (
            ast.literal_eval(node.value)
            for node in ast.walk(module)
            if isinstance(node, ast.Assign) and node.targets[0].id == x
        )

        try:
            return (
                next(itr("__author__")),
                next(itr("__email__")),
                next(itr("__license__")),
            )
        except StopIteration as stop_exception:
            raise ValueError(
                "One of author, email, license cannot be found in {}".format(filename)
            ) from stop_exception

    def missing_requirements(self, specifiers):
        """Find what's missing"""
        for specifier in specifiers:
            try:
                pkg_resources.require(specifier)
            except pkg_resources.DistributionNotFound:
                yield specifier

    def install_requirements(self, requires):
        """Install the listed requirements"""
        # Temporarily install dependencies required by setup.py before trying to import them.
        sys.path[0:0] = ["setup-requires"]
        pkg_resources.working_set.add_entry("setup-requires")

        to_install = list(self.missing_requirements(requires))
        if to_install:
            cmd = [
                sys.executable,
                "-m",
                "pip",
                "install",
                "-t",
                "setup-requires",
            ] + to_install
            subprocess.call(cmd)

    def get_long_description(self, filename="README.md"):
        """I really prefer Markdown to reStructuredText. PyPi does not."""
        return open(filename).read()
