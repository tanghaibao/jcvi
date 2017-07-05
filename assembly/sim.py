#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Simulate Illumina sequencing reads.
"""

import sys
import logging
import math

from jcvi.formats.fasta import Fasta
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


def main():

    actions = (
        ('wgsim', 'sample paired end reads using dwgsim'),
        ('eagle', 'simulate Illumina reads using EAGLE'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def eagle(args):
    """
    %prog eagle fastafile

    """
    p = OptionParser(eagle.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())


def wgsim(args):
    """
    %prog wgsim fastafile

    Run dwgsim on fastafile.
    """
    p = OptionParser(wgsim.__doc__)
    p.add_option("--erate", default=.01, type="float",
                 help="Base error rate of the read [default: %default]")
    p.add_option("--noerrors", default=False, action="store_true",
                 help="Simulate reads with no errors [default: %default]")
    p.add_option("--distance", default=500, type="int",
                 help="Outer distance between the two ends [default: %default]")
    p.add_option("--genomesize", type="int",
                 help="Genome size in Mb [default: estimate from data]")
    p.add_option("--readlen", default=150, type="int",
                 help="Length of the read")
    p.set_outfile(outfile=None)
    p.set_depth(depth=10)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    pf = fastafile.split(".")[0]

    genomesize = opts.genomesize
    size = genomesize * 1000000 if genomesize else Fasta(fastafile).totalsize
    depth = opts.depth
    readlen = opts.readlen
    readnum = int(math.ceil(size * depth / (2 * readlen)))

    distance = opts.distance
    stdev = distance / 10

    outpf = opts.outfile or "{0}.{1}bp.{2}x".format(pf, distance, depth)

    logging.debug("Total genome size: {0} bp".format(size))
    logging.debug("Target depth: {0}x".format(depth))
    logging.debug("Number of read pairs (2x{0}): {1}".format(readlen, readnum))

    if opts.noerrors:
        opts.erate = 0

    cmd = "dwgsim -e {0} -E {0}".format(opts.erate)
    if opts.noerrors:
        cmd += " -r 0 -R 0 -X 0 -y 0"

    cmd += " -d {0} -s {1}".format(distance, stdev)
    cmd += " -N {0} -1 {1} -2 {1}".format(readnum, readlen)
    cmd += " {0} {1}".format(fastafile, outpf)
    sh(cmd)


if __name__ == '__main__':
    main()
