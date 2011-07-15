#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Using CD-HIT (CD-HIT-454) in particular to remove duplicate reads.
"""

import sys

from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.base import read_block
from jcvi.apps.base import ActionDispatcher, set_grid, debug, sh
debug()


CDPATH="~/htang/export/cd-hit-v4.5.5-2011-03-31/"


def main():

    actions = (
        ('deduplicate', 'use `cd-hit-454` to remove duplicate reads'),
        ('summary', 'parse cdhit.clstr file to get distribution of cluster sizes'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def summary(args):
    """
    %prog summary cdhit.clstr

    Parse cdhit.clstr file to get distribution of cluster sizes.
    """
    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    from math import log

    clstrfile, = args
    bins = defaultdict(int)
    fp = open(clstrfile)
    unique = 0
    total = 0
    for clstr, members in read_block(fp, ">"):
        size = len(members)
        unique += 1
        total += size
        log2size = int(log(size, 2))
        bins[log2size] += 1

    # Print out a distribution
    print >> sys.stderr, "Total: {0}, Unique: {1}".format(total, unique)
    for size, number in sorted(bins.items()):
        lb, ub = 2 ** size, 2 ** (size + 1)
        print >> sys.stderr, "Counts in range [{0}, {1}): {2}".format(lb, ub, number)


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


if __name__ == '__main__':
    main()
