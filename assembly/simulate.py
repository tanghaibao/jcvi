#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Simulate fake reads from genome for benchmarking.
"""

import sys
import logging

from optparse import OptionParser

from jcvi.formats.fasta import Fasta
from jcvi.apps.base import ActionDispatcher, debug, sh
debug()


def main():

    actions = (
        ('wgsim', 'sample paired end reads using dwgsim'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def wgsim(args):
    """
    %prog wgsim fastafile outprefix

    Run dwgsim on fastafile.
    """
    p = OptionParser(wgsim.__doc__)
    p.add_option("--erate", default=.02, type="float",
                 help="Base error rate of the read [default: %default]")
    p.add_option("--distance", default=500, type="int",
                 help="Outer distance between the two ends [default: %default]")
    p.add_option("--depth", default=50, type="int",
                 help="Target depth (aka base coverage) [default: %default]")
    p.add_option("--readlen", default=100, type="int",
                 help="Length of the read [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, outpf = args
    size = Fasta(fastafile).totalsize
    depth = opts.depth
    readlen = opts.readlen
    readnum = size * depth / (2 * readlen)

    distance = opts.distance
    stdev = distance / 5
    distance -= 2 * readlen  # Outer distance => Inner distance

    logging.debug("Total FASTA size: {0} bp".format(size))
    logging.debug("Target depth: {0}x".format(depth))
    logging.debug("Number of read pairs (2x{0}): {1}".format(readlen, readnum))

    cmd = "dwgsim -e {0} -E {0}".format(opts.erate)
    cmd += " -d {0} -s {1}".format(distance, stdev)
    cmd += " -N {0} -1 {1} -2 {1}".format(readnum, readlen)
    cmd += " {0} {1}".format(fastafile, outpf)
    sh(cmd)


if __name__ == '__main__':
    main()
