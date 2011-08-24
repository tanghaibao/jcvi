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

from jcvi.formats.agp import AGP
from jcvi.formats.fasta import Fasta, SeqIO
from jcvi.formats.blast import BlastLine
from jcvi.formats.coords import Overlap_types
from jcvi.apps.entrez import fetch
from jcvi.apps.base import ActionDispatcher, debug, popen, mkdir, sh
debug()


class Overlap (object):

    def __init__(self, aid, bid, otype, pctid, hitlen, orientation):
        self.aid = aid.split("|")[3] if aid.count("|") >= 3 else aid
        self.bid = bid.split("|")[3] if bid.count("|") >= 3 else bid
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
        ('bes', 'confirm the BES mapping'),
        ('flip', 'flip the FASTA sequences according to a set of references'),
        ('overlap', 'check terminal overlaps between two records'),
        ('neighbor', 'check neighbors of a component in agpfile'),
        ('add', 'add a component into agpfile'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def add(args):
    """
    %prog add allfasta clonename

    Insert a component into agpfile by aligning to the best hit in pool and see
    if they have good overlaps.
    """
    from jcvi.apps.command import run_megablast

    p = OptionParser(add.__doc__)
    p.add_option("-n", type="int", default=1,
            help="Take best N hits [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    allfasta, clonename = args
    fastadir = "fasta"
    fetch([clonename, "--skipcheck", "--outdir=" + fastadir])

    infile = op.join(fastadir, clonename + ".fasta")
    outfile = "{0}.{1}.blast".format(clonename, allfasta.split(".")[0])
    run_megablast(infile=infile, outfile=outfile, db=allfasta, \
            pctid=99, hitlen=1000)

    blasts = [BlastLine(x) for x in open(outfile)]
    besthit = set()
    for b in blasts:
        if b.query.count("|") >= 3:
            b.query = b.query.split("|")[3]

        if b.subject.count("|") >= 3:
            b.subject = b.subject.split("|")[3]

        b.query = b.query.rsplit(".", 1)[0]
        b.subject = b.subject.rsplit(".", 1)[0]

        if b.query == b.subject:
            continue

        besthit.add(b.subject)
        if len(besthit) == opts.n:
            break

    for b in besthit:
        overlap([clonename, b, "--dir=" + fastadir])


def bes(args):
    """
    %prog bes bacfasta clonename

    Use the clone name to download BES gss sequences from Genbank, map and then
    visualize.
    """
    from jcvi.apps.command import run_blat

    p = OptionParser(bes.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bacfasta, clonename = args

    fetch([clonename, "--database=nucgss", "--skipcheck"])
    besfasta = clonename + ".fasta"
    blatfile = clonename + ".bes.blat"
    run_blat(infile=besfasta, outfile=blatfile, db=bacfasta, \
             pctid=95, hitlen=100)

    aid, asize = Fasta(bacfasta).itersizes().next()

    width = 50
    msg = "=" * width
    msg += "  " + aid
    print >> sys.stderr, msg

    ratio = width * 1. / asize
    _ = lambda x: int(round(x * ratio, 0))
    blasts = [BlastLine(x) for x in open(blatfile)]
    for b in blasts:
        if b.orientation == '+':
            msg = " " * _(b.sstart) + "->"
        else:
            msg = " " * (_(b.sstop) - 2) + "<-"
        msg += " " * (width - len(msg) + 2)
        msg += b.query
        if b.orientation == '+':
            msg += " (hang={0})".format(b.sstart - 1)
        else:
            msg += " (hang={0})".format(asize - b.sstop)

        print >> sys.stderr, msg


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
    from jcvi.apps.command import BLPATH

    p = OptionParser(overlap.__doc__)
    p.add_option("--dir", default=os.getcwd(),
            help="Download sequences to dir [default: %default]")
    p.add_option("--qreverse", default=False, action="store_true",
            help="Reverse seq a [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    afasta, bfasta = args
    dir = opts.dir

    if not op.exists(afasta):
        fetch([afasta, "--skipcheck", "--outdir=" + dir])
        afasta += ".fasta"

    if not op.exists(bfasta):
        fetch([bfasta, "--skipcheck", "--outdir=" + dir])
        bfasta += ".fasta"

    afasta = op.join(dir, afasta)
    bfasta = op.join(dir, bfasta)

    cmd = BLPATH("blastn")
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
    otype = b.overlap(asize, bsize, qreverse=opts.qreverse)
    o = Overlap(aid, bid, otype, pctid, hitlen, orientation)
    print >> sys.stderr, str(o)

    return o


def neighbor(args):
    """
    %prog neighbor agpfile componentID

    Check overlaps of a particular component in agpfile.
    """
    p = OptionParser(neighbor.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    agpfile, componentID = args
    fastadir = "fasta"

    cmd = "grep"
    cmd += " --color -C2 {0} {1}".format(componentID, agpfile)
    sh(cmd)

    agp = AGP(agpfile)
    aorder = agp.order
    i, c = aorder[componentID]
    north, south = agp.getNorthSouthClone(i)

    if not north.isCloneGap:
        ar = [north.component_id, componentID, "--dir=" + fastadir]
        if north.orientation == '-':
            ar += ["--qreverse"]
        overlap(ar)

    if not south.isCloneGap:
        ar = [componentID, south.component_id, "--dir=" + fastadir]
        if c.orientation == '-':
            ar += ["--qreverse"]
        overlap(ar)


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


if __name__ == '__main__':
    main()
