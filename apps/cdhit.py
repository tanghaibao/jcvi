#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Using CD-HIT (CD-HIT-454) in particular to remove duplicate reads.
"""

import sys
import logging

from math import log
from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.base import read_block
from jcvi.utils.cbook import percentage
from jcvi.apps.base import ActionDispatcher, set_grid, debug, sh, CDPATH
debug()


def main():

    actions = (
        ('ids', 'get the representative ids from clstr file'),
        ('deduplicate', 'use `cd-hit-454` to remove duplicate reads'),
        ('summary', 'parse cdhit.clstr file to get distribution of cluster sizes'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def ids(args):
    """
    %prog ids cdhit.clstr

    Get the representative ids from clstr file.
    """
    p = OptionParser(ids.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clstrfile, = args
    assert clstrfile.endswith(".clstr")

    idsfile = clstrfile.replace(".clstr", ".ids")
    fp = open(clstrfile)
    fw = open(idsfile, "w")
    nreads = 0
    for row in fp:
        if row[0] == '>':
            continue
        name = row.split('>', 1)[1].split("...")[0]
        if row.rstrip()[-1] == '*':
            print >> fw, name
            nreads += 1

    logging.debug("A total of {0} unique reads written to `{1}`.".\
            format(nreads, idsfile))


def summary(args):
    """
    %prog summary cdhit.clstr

    Parse cdhit.clstr file to get distribution of cluster sizes.
    """
    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clstrfile, = args
    assert clstrfile.endswith(".clstr")

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

    from jcvi.graphics.histogram import loghistogram
    # Print out a distribution
    print >> sys.stderr, "Unique: {0}".format(percentage(unique, total))
    loghistogram(bins)


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
