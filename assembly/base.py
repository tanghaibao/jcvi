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
from jcvi.apps.base import ActionDispatcher, debug, CAPATH
debug()


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
    """
    p = OptionParser(n50.__doc__)

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(p.print_help())

    ctgsizes = []
    bins = defaultdict(int)
    for filename in args:
        for row in open(filename):
            try:
                ctgsize = int(row.strip())
            except ValueError:
                continue
            log2ctgsize = int(log(ctgsize, 2))
            bins[log2ctgsize] += 1
            ctgsizes.append(ctgsize)

    a50, l50, nn50 = calculate_A50(ctgsizes)
    sumsize = sum(ctgsizes)
    maxsize = max(ctgsizes)
    print >> sys.stderr, ", ".join(args)
    print >> sys.stderr, "Length={1} L50={2} Max={3}".\
            format(filename, sumsize, l50, maxsize)
    loghistogram(bins)


def main():

    actions = (
            ('n50', "Given a list of numbers, calculate N50"),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == '__main__':
    main()
