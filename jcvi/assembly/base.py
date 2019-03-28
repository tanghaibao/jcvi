#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Base utilties for genome assembly related calculations and manipulations
"""
from __future__ import print_function

import os.path as op
import sys

from math import log
ln2 = log(2)

import numpy as np
from bisect import bisect

from jcvi.formats.base import must_open
from jcvi.formats.fasta import Fasta
from jcvi.apps.base import OptionParser, ActionDispatcher, glob


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
        self.type = types[pf]
        self.stddev = size / 6 if self.type == "jumping" else size / 9
        self.paired = 0 if size == 0 else 1
        self.read_orientation = "outward" if pf == "MP" else "inward"
        self.reverse_seq = 1 if pf == "MP" else 0
        self.asm_flags = 3 if pf != "MP" else 2
        if not self.paired:
            self.read_orientation = ""

    def get_lib_seq(self, wildcard, prefix, readlen, rank):
        # lib_seq wildcard prefix insAvg insSdev avgReadLen hasInnieArtifact
        # isRevComped useForContigging scaffRound useForGapClosing 5pWiggleRoom
        # 3pWiggleRoom (used by MERACULOUS)
        useForContigging = useForGapClosing = int(self.asm_flags == 3)
        return ("lib_seq", wildcard, prefix, self.size,
                self.stddev, readlen, int(self.type == "jumping"),
                self.reverse_seq, useForContigging, rank,
                useForGapClosing, 0, 0)


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
                ctgsize = int(float(row.split()[-1]))
            except ValueError:
                continue
            ctgsizes.append(ctgsize)

    a50, l50, nn50 = calculate_A50(ctgsizes)
    sumsize = sum(ctgsizes)
    minsize = min(ctgsizes)
    maxsize = max(ctgsizes)
    n = len(ctgsizes)
    print(", ".join(args), file=sys.stderr)

    summary = (sumsize, l50, nn50, minsize, maxsize, n)
    print(" ".join("{0}={1}".format(a, b) for a, b in \
                        zip(header, summary)), file=sys.stderr)
    loghistogram(ctgsizes)

    if opts.print0:
        print("\t".join(str(x) for x in (",".join(args), sumsize, l50)))

    return zip(header, summary)


if __name__ == '__main__':
    main()
