#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Locate IES sequences within MIC genome of tetrahymena.
"""

import sys

from jcvi.formats.bed import Bed, sort
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update


def main():

    actions = (
        ('deletion', 'find IES based on mapping MAC reads'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def deletion(args):
    """
    %prog deletion mac.mic.bed

    Find IES based on mapping MAC reads to MIC genome.
    """
    p = OptionParser(deletion.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    sortedbedfile = bedfile.rsplit(".", 1)[0] + ".sorted.bed"
    if need_update(bedfile, sortedbedfile):
        sort([bedfile, "--accn"])

    bed = Bed(sortedbedfile)


if __name__ == '__main__':
    main()
