#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Base utilties for genome assembly related calculations and manipulations
"""

import os.path as op
import sys
import logging

from math import log
ln2 = log(2)

import numpy as np
from bisect import bisect

from jcvi.formats.base import must_open
from jcvi.formats.fasta import Fasta
from jcvi.apps.base import OptionParser, ActionDispatcher, glob, sh


types = {"PE": "fragment", "MP": "jumping", "TT": "jumping", "LL": "long"}
header = ("Length", "L50", "N50", "Min", "Max", "N")

FastqNamings = """
    The naming schemes for the fastq files are.

    PE-376.fastq (paired end)
    MP-3000.fastq (mate pairs)
    TT-3000.fastq (mate pairs, but from 454 data, so expected to be +-)
    LL-0.fastq (long reads)

    Paired reads in different files must be in the form of (note the .1. and .2.):
    PE-376.1.fastq and PE-376.2.fastq to be considered

    The reads are assumed to be NOT paired if the number after the PE-, MP-,
    etc. is 0. Otherwise, they are considered paired at the given distance.
"""


class Library (object):
    """
    The sequence files define a library.
    """
    def __init__(self, library_name):

        self.library_name = library_name
        if "-" in library_name:
            pf, size = library_name.split("-", 1)
            assert pf in types, \
                "Library prefix must be one of {0}".format(types.keys())
        else:
            pf, size = "PE", 0

        self.size = size = int(size)
        self.stddev = size / 5
        self.type = types[pf]
        self.paired = 0 if size == 0 else 1
        self.read_orientation = "outward" if pf == "MP" else "inward"
        self.reverse_seq = 1 if pf == "MP" else 0
        self.asm_flags = 3 if pf != "MP" else 2
        if not self.paired:
            self.read_orientation = ""


def get_libs(args):
    from itertools import groupby

    fnames = args or glob("*.fastq*")
    fnames = sorted(fnames)
    for x in fnames:
        assert op.exists(x), "File `{0}` not found.".format(x)

    library_name = lambda x: "-".join(\
                op.basename(x).split(".")[0].split("-")[:2])
    libs = [(Library(x), sorted(fs)) for x, fs in \
                groupby(fnames, key=library_name)]

    libs.sort(key=lambda x: x[0].size)
    return libs


def calculate_A50(ctgsizes, cutoff=0, percent=50):
    """
    Given an array of contig sizes, produce A50, N50, and L50 values
    """

    ctgsizes = np.array(ctgsizes, dtype="int")
    ctgsizes = np.sort(ctgsizes)[::-1]
    ctgsizes = ctgsizes[ctgsizes >= cutoff]

    a50 = np.cumsum(ctgsizes)

    total = np.sum(ctgsizes)
    idx = bisect(a50, total * percent / 100.)
    l50 = ctgsizes[idx]
    n50 = idx + 1

    return a50, l50, n50


"""
Discriminator A-statistics:

If n reads are uniform sample of the genome of length G,
we expect k = n * delta / G to start in a region of length delta

Use poisson distribution:
A(delta, k) = ln(prob(1-copy) / prob(2-copies)) = n * delta / G - k * ln2
"""


def Astat(delta, k, G, n):
    """
    delta: contig size
    k: reads mapped in contig
    G: total genome size
    n: total reads mapped to genome
    """
    return n * delta * 1. / G - k * ln2


def main():

    actions = (
            ('n50', "Given FASTA or a list of contig sizes, calculate N50"),
            ('wgsim', 'sample paired end reads using dwgsim'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def n50(args):
    """
    %prog n50 filename

    Given a file with a list of numbers denoting contig lengths, calculate N50.
    Input file can be both FASTA or a list of sizes.
    """
    from jcvi.graphics.histogram import loghistogram

    p = OptionParser(n50.__doc__)
    p.add_option("--print0", default=False, action="store_true",
                 help="Print size and L50 to stdout [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    ctgsizes = []

    # Guess file format
    probe = open(args[0]).readline()[0]
    isFasta = (probe == '>')
    if isFasta:
        for filename in args:
            f = Fasta(filename)
            ctgsizes += list(b for a, b in f.itersizes())

    else:
        for row in must_open(args):
            try:
                ctgsize = int(row.split()[-1])
            except ValueError:
                continue
            ctgsizes.append(ctgsize)

    a50, l50, nn50 = calculate_A50(ctgsizes)
    sumsize = sum(ctgsizes)
    minsize = min(ctgsizes)
    maxsize = max(ctgsizes)
    n = len(ctgsizes)
    print >> sys.stderr, ", ".join(args)

    summary = (sumsize, l50, nn50, minsize, maxsize, n)
    print >> sys.stderr, " ".join("{0}={1}".format(a, b) for a, b in \
                        zip(header, summary))
    loghistogram(ctgsizes)

    if opts.print0:
        print "\t".join(str(x) for x in (",".join(args), sumsize, l50))

    return zip(header, summary)


def wgsim(args):
    """
    %prog wgsim fastafile

    Run dwgsim on fastafile.
    """
    p = OptionParser(wgsim.__doc__)
    p.add_option("--erate", default=.02, type="float",
                 help="Base error rate of the read [default: %default]")
    p.add_option("--distance", default=500, type="int",
                 help="Outer distance between the two ends [default: %default]")
    p.add_option("--genomesize", type="int",
                 help="Genome size in Mb [default: estimate from data]")
    p.add_option("--readlen", default=100, type="int",
                 help="Length of the read [default: %default]")
    p.add_option("--noerrors", default=False, action="store_true",
                 help="Simulate reads with no errors [default: %default]")
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
    readnum = size * depth / (2 * readlen)

    distance = opts.distance
    stdev = distance / 5

    outpf = "{0}.{1}bp.{2}x".format(pf, distance, depth)
    distance -= 2 * readlen  # Outer distance => Inner distance
    assert distance >= 0, "Outer distance must be >= 2 * readlen"

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
