#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Optical map alignment parser.
"""

import sys

from optparse import OptionParser
from xml.etree.ElementTree import ElementTree

from jcvi.apps.base import ActionDispatcher, debug
debug()


class OpticalMap (object):

    def __init__(self, xmlfile):
        tree = ElementTree()
        self.root = tree.parse(xmlfile)

    def iter_maps(self):
        for e in self.root.findall("restriction_map"):
            yield e

    def iter_alignments(self):
        for e in self.root.findall("map_alignment"):
            yield e


def main():

    actions = (
        ('summary', 'print summary of optical map alignments'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def summary(args):
    """
    %prog summary xmlfile

    Print summary of optical map alignment.
    """
    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    xmlfile, = args
    om = OpticalMap(xmlfile)
    print list(om.iter_maps())
    print list(om.iter_alignments())


if __name__ == '__main__':
    main()
