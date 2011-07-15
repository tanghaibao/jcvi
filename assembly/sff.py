#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Deal with Roche SFF file format, including converting to FASTA, split linkers,
and remove duplicates.
"""

import sys

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, set_grid, debug, sh
debug()


CDPATH="~/htang/export/cd-hit-v4.5.5-2011-03-31/"


def main():

    actions = (
        ('fasta', 'use `sffinfo` to convert to fasta format'),
        ('deduplicate', 'use `cd-hit-454` to remove duplicate reads'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def deduplicate(args):
    """
    %prog deduplicate fastafile

    Wraps `cd-hit-454` to remove duplicate reads.
    """
    p = OptionParser(deduplicate.__doc__)
    set_grid(p)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args

    cmd = CDPATH + "cd-hit-454 -M 0 -T 0 " + \
            "-i {0} -o {0}.cdhit".format(fastafile)
    sh(cmd, grid=opts.grid)


def fasta(args):
    """
    %prog fasta sffiles

    Wraps `sffinfo` to convert sffile to fastafile, can run on grid.
    """
    p = OptionParser(fasta.__doc__)
    set_grid(p)

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    for sffile in args:
        assert sffile.endswith(".sff")

        fastafile = sffile.replace(".sff", ".fasta")
        cmd = "sffinfo -seq {0}".format(sffile)

        sh(cmd, outfile=fastafile, grid=opts.grid)


if __name__ == '__main__':
    main()
