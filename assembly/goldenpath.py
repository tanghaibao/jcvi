#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Procedures to validate and update golden path of a genome assembly. This relies
heavily on formats.agp, and further includes several algorithms, e.g. overlap
detection.
"""

import os
import os.path as op
import sys

from optparse import OptionParser

from jcvi.formats.fasta import Fasta, SeqIO
from jcvi.formats.blast import BlastLine
from jcvi.formats.coords import Overlap_types
from jcvi.apps.entrez import fetch
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
        ('flip', 'flip the FASTA sequences according to a set of references'),
        ('overlap', 'check terminal overlaps between two records'),
        ('overlapbatch', 'check overlaps between adjacent components in agpfile'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def flip(args):
    """
    %prog flip fastafile

    Go through each FASTA record, check against Genbank file and determines
    whether or not to flip the sequence. This is useful before updates of the
    sequences to make sure the same orientation is used.
    """
    p = OptionParser(flip.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    outfastafile = fastafile.rsplit(".", 1)[0] + ".flipped.fasta"
    fo = open(outfastafile, "w")
    f = Fasta(fastafile, lazy=True)
    for name, rec in f.iteritems_ordered():
        tmpfasta = "a.fasta"
        fw = open(tmpfasta, "w")
        SeqIO.write([rec], fw, "fasta")
        fw.close()

        o = overlap([tmpfasta, name])
        if o.orientation == '-':
            rec.seq = rec.seq.reverse_complement()

        SeqIO.write([rec], fo, "fasta")
        os.remove(tmpfasta)


def overlap(args):
    """
    %prog overlap <a|a.fasta> <b|b.fasta>

    Check overlaps between two fasta records. The arguments can be genBank IDs
    instead of FASTA files. In case of IDs, the sequences will be downloaded
    first.
    """
    p = OptionParser(overlap.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    afasta, bfasta = args

    if not op.exists(afasta):
        fetch([afasta])
        afasta += ".fasta"

    if not op.exists(bfasta):
        fetch([bfasta])
        bfasta += ".fasta"

    cmd = "blastn"
    cmd += " -query {0} -subject {1}".format(afasta, bfasta)
    cmd += " -evalue 0.01 -outfmt 6"

    fp = popen(cmd)
    hsps = fp.readlines()
    if len(hsps) == 0:
        print >> sys.stderr, "No match found."
        return None

    besthsp = hsps[0]
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
