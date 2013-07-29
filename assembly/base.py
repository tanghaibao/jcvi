#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Base utilties for genome assembly related calculations and manipulations
"""

import os.path as op
import sys
from math import log
ln2 = log(2)

import numpy as np
from bisect import bisect
from optparse import OptionParser

from jcvi.graphics.histogram import loghistogram
from jcvi.formats.base import must_open
from jcvi.formats.fasta import Fasta
from jcvi.apps.base import ActionDispatcher, debug
debug()


orientationlabels = {"++": "normal", "+-": "innie", "-+": "outie", "--": "antinormal"}
orientationflips = {"++": "--", "+-": "-+", "-+": "+-", "--": "++"}
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
    loghistogram(ctgsizes, summary=False)

    if opts.print0:
        print "\t".join(str(x) for x in (",".join(args), sumsize, l50))

    return zip(header, summary)


def main():

    actions = (
            ('n50', "Given FASTA or a list of contig sizes, calculate N50"),
            ('allstats', "Summarize multiple FASTA in a table"),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def allstats(args):
    """
    %prog allstats fastafiles

    Summarize multiple FASTA in a table.
    """
    from jcvi.utils.table import tabulate

    p = OptionParser(allstats.__doc__)
    p.add_option("--exclude", help="Exclude statistics, must be {0}, "
                      "multiple separated by comma [default: %default]".\
                      format("|".join(header))
                 )

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fastafiles = args
    exclude = opts.exclude.split(",")
    assert all(x in header for x in exclude)

    tabledict = {}
    for fastafile in fastafiles:
        pf = fastafile.rsplit(".", 1)[0]
        for key, val in n50([fastafile]):
            if key in exclude:
                continue
            tabledict[(pf, key)] = val

    table = tabulate(tabledict)
    print >> sys.stderr, table


if __name__ == '__main__':
    main()
