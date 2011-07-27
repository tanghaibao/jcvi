#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Base utilties for genome assembly related calculations and manipulations
"""

import sys
from math import log
ln2 = log(2)

import numpy as np
from bisect import bisect
from collections import defaultdict
from optparse import OptionParser

from jcvi.graphics.histogram import loghistogram
from jcvi.formats.base import must_open
from jcvi.formats.fasta import Fasta
from jcvi.apps.command import CAPATH
from jcvi.apps.base import ActionDispatcher, debug
debug()


orientationlabels = {"++": "normal", "+-": "innie", "-+": "outie", "--": "antinormal"}


def calculate_A50(ctgsizes, cutoff=0):
    """
    Given an array of contig sizes, produce A50, N50, and L50 values
    """

    ctgsizes = np.array(ctgsizes, dtype="int")
    ctgsizes = np.sort(ctgsizes)[::-1]
    ctgsizes = ctgsizes[ctgsizes >= cutoff]

    a50 = np.cumsum(ctgsizes)

    total = np.sum(ctgsizes)
    idx = bisect(a50, total / 2)
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


def n50(args):
    """
    %prog n50 filename

    Given a file with a list of numbers denoting contig lengths, calculate N50.
    Input file can be both FASTA or a list of sizes.
    """
    p = OptionParser(n50.__doc__)

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    ctgsizes = []
    bins = defaultdict(int)

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
                ctgsize = int(row.strip())
            except ValueError:
                continue
            ctgsizes.append(ctgsize)

    for ctgsize in ctgsizes:
        log2ctgsize = int(log(ctgsize, 2))
        bins[log2ctgsize] += 1

    a50, l50, nn50 = calculate_A50(ctgsizes)
    sumsize = sum(ctgsizes)
    minsize = min(ctgsizes)
    maxsize = max(ctgsizes)
    print >> sys.stderr, ", ".join(args)
    print >> sys.stderr, "Length={0} L50={1} Min={2} Max={3} N={4}".\
            format(sumsize, l50, minsize, maxsize, len(ctgsizes))
    loghistogram(bins)


def main():

    actions = (
            ('n50', "Given a list of numbers, calculate N50"),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == '__main__':
    main()
