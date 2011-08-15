#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Procedures to validate and update golden path of a genome assembly. This relies
heavily on formats.agp, and further includes several algorithms, e.g. overlap
detection.
"""

import sys

from optparse import OptionParser

from jcvi.formats.fasta import Fasta
from jcvi.formats.blast import BlastLine
from jcvi.formats.coords import Overlap_types
from jcvi.apps.base import ActionDispatcher, debug, popen, mkdir
debug()


class Overlap (object):

    def __init__(self, aid, bid, otype, pctid, hitlen, orientation):
        self.aid = aid
        self.bid = bid
        self.otype = otype

        self.pctid = pctid
        self.hitlen = hitlen
        self.orientation = orientation

    def __str__(self):
        ov = Overlap_types[self.otype]
        s = "{0} - {1}: {2} ".format(self.aid, self.bid, ov)
        s += "Overlap: {0} Identity: {1}% Orientation: {2}".\
            format(self.hitlen, self.pctid, self.orientation)
        return s

    def isTerminal(self, length_cutoff=2000, pctid_cutoff=99):
        return self.otype in (1, 2)

    def isGoodQuality(self, length_cutoff=2000, pctid_cutoff=99):
        return self.hitlen >= length_cutoff and \
               self.pctid >= pctid_cutoff


def main():

    actions = (
        ('overlap', 'check terminal overlaps between two records'),
        ('overlapbatch', 'check overlaps between adjacent components in agpfile'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def overlap(args):
    """
    %prog overlap a.fasta b.fasta

    Check overlaps between two fasta records.
    """
    p = OptionParser(overlap.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    afasta, bfasta = args

    cmd = "blastn"
    cmd += " -query {0} -subject {1}".format(afasta, bfasta)
    cmd += " -evalue 0.01 -outfmt 6"

    fp = popen(cmd)
    besthsp = fp.readline()
    b = BlastLine(besthsp)
    pctid = b.pctid
    hitlen = b.hitlen
    orientation = b.orientation

    aid, asize = Fasta(afasta).itersizes().next()
    bid, bsize = Fasta(bfasta).itersizes().next()
    otype = b.overlap(asize, bsize)
    o = Overlap(aid, bid, otype, pctid, hitlen, orientation)
    print >> sys.stderr, str(o)

    return o


def overlapbatch(args):
    """
    %prog overlapbatch agpfile componentfasta

    Check overlaps between adjacent components in an AGP file.
    """
    p = OptionParser(overlapbatch.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    agpfile, componentfasta = args

    fastadir = "fasta"
    mkdir(fastadir, overwrite=False)

    cmd = "faSplit"
    cmd += " byname {0} {1}/".format(componentfasta, fastadir)
    sh(cmd)


if __name__ == '__main__':
    main()
