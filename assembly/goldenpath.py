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
import shutil
import logging

from optparse import OptionParser
from itertools import groupby

from jcvi.formats.agp import AGP, TPF, get_phase
from jcvi.formats.base import BaseFile
from jcvi.formats.fasta import Fasta, SeqIO
from jcvi.formats.blast import BlastLine
from jcvi.formats.coords import Overlap_types
from jcvi.formats.base import must_open
from jcvi.utils.cbook import memoized
from jcvi.apps.entrez import fetch
from jcvi.apps.base import ActionDispatcher, debug, popen, mkdir, sh
debug()


class Overlap (object):

    def __init__(self, blastline, asize, bsize):
        b = blastline
        aid = b.query
        bid = b.subject

        self.aid = aid.split("|")[3] if aid.count("|") >= 3 else aid
        self.bid = bid.split("|")[3] if bid.count("|") >= 3 else bid
        self.asize = asize
        self.bsize = bsize

        self.qstart = b.qstart
        self.qstop = b.qstop
        self.sstart = b.sstart
        self.sstop = b.sstop

        self.pctid = b.pctid
        self.hitlen = b.hitlen
        self.orientation = b.orientation

        self.otype = self.get_otype()

    def __str__(self):
        ov = Overlap_types[self.otype]
        s = "{0} - {1}: {2} ".format(self.aid, self.bid, ov)
        s += "Overlap: {0} Identity: {1}% Orientation: {2}".\
            format(self.hitlen, self.pctid, self.orientation)
        return s

    @property
    def certificateline(self):
        terminal_tag = "Terminal" if self.isTerminal() else "Non-terminal"
        return "\t".join(str(x) for x in (self.bid, \
                          self.asize, self.qstart, self.qstop, \
                          self.orientation, terminal_tag))

    def isTerminal(self, length_cutoff=100, pctid_cutoff=99):
        return self.isGoodQuality(length_cutoff, pctid_cutoff) \
               and self.otype in (1, 2)

    def isGoodQuality(self, length_cutoff=100, pctid_cutoff=99):
        return self.hitlen >= length_cutoff and \
               self.pctid >= pctid_cutoff

    def get_hangs(self):
        """
        Determine the type of overlap given query, ref alignment coordinates
        Consider the following alignment between sequence a and b:

        aLhang \              / aRhang
                \------------/
                /------------\
        bLhang /              \ bRhang

        Terminal overlap: a before b, b before a
        Contain overlap: a in b, b in a
        """
        aLhang, aRhang = self.qstart - 1, self.asize - self.qstop
        bLhang, bRhang = self.sstart - 1, self.bsize - self.sstop
        if self.orientation == '-':
            bLhang, bRhang = bRhang, bLhang

        return aLhang, aRhang, bLhang, bRhang

    def print_graphic(self, qreverse=False):
        """
        >>>>>>>>>>>>>>>>>>>             seqA (alen)
                  ||||||||
                 <<<<<<<<<<<<<<<<<<<<<  seqB (blen)
        """
        aLhang, aRhang, bLhang, bRhang = self.get_hangs()

        if qreverse:
            aLhang, aRhang = aRhang, aLhang
            bLhang, bRhang = bRhang, bLhang

        achar = ">"
        bchar = "<" if self.orientation == '-' else ">"
        if qreverse:
            achar = "<"
            bchar = {">" : "<", "<" : ">"}[bchar]

        print >> sys.stderr, aLhang, aRhang, bLhang, bRhang
        width = 50  # Canvas
        hitlen = self.hitlen
        lmax = max(aLhang, bLhang)
        rmax = max(aRhang, bRhang)
        bpwidth = lmax + hitlen + rmax
        ratio = width * 1. / bpwidth

        _ = lambda x: int(round(x * ratio, 0))
        a1, a2 = _(aLhang), _(aRhang)
        b1, b2 = _(bLhang), _(bRhang)
        hit = max(_(hitlen), 1)

        msg = " " * max(b1 - a1, 0)
        msg += achar * (a1 + hit + a2)
        msg += " " * (width - len(msg) + 2)
        msg += "{0} ({1})".format(self.aid, self.asize)
        print >> sys.stderr, msg

        msg = " " * max(a1, b1)
        msg += "|" * hit
        print >> sys.stderr, msg

        msg = " " * max(a1 - b1, 0)
        msg += bchar * (b1 + hit + b2)
        msg += " " * (width - len(msg) + 2)
        msg += "{0} ({1})".format(self.bid, self.bsize)
        print >> sys.stderr, msg

    def get_otype(self, max_hang=10000):
        aLhang, aRhang, bLhang, bRhang = self.get_hangs()

        s1 = aLhang + bRhang
        s2 = aRhang + bLhang
        s3 = aLhang + aRhang
        s4 = bLhang + bRhang

        # Dovetail (terminal) overlap
        if s1 < max_hang:
            type = 2  # b ~ a
        elif s2 < max_hang:
            type = 1  # a ~ b
        # Containment overlap
        elif s3 < max_hang:
            type = 3  # a in b
        elif s4 < max_hang:
            type = 4  # b in a
        else:
            type = 0

        return type


class CertificateLine (object):

    """
    North  chr1  2  0  AC229737.8  telomere     58443
    South  chr1  2  1  AC229737.8  AC202463.29  58443  37835  58443  + Non-terminal
    """

    def __init__(self, line):
        args = line.split()
        self.tag = args[0]
        self.chr = args[1]
        self.aphase = int(args[2])
        self.bphase = int(args[3])
        self.aid = args[4]
        self.bid = args[5]
        self.asize = int(args[6])

        if len(args) == 7:
            self.is_gap = True
            return

        self.is_gap = False

        if len(args) == 8:
            assert args[7] == "None"
            self.terminal = "Non-terminal"
            return

        self.astart = int(args[7])
        self.astop = int(args[8])
        self.orientation = args[9]
        self.terminal = args[10]

    def isTerminal(self):
        return self.terminal == "Terminal"

    def __str__(self):
        ar = [self.tag, self.chr, self.aphase, self.bphase, self.aid, self.bid]
        if not self.is_gap:

            ar += [self.asize, self.astart, self.astop,
                   self.orientation, self.terminal]

        return "\t".join(str(x) for x in ar)


class Certificate (BaseFile):

    gapsize = 100000
    gaps = dict(telomere=gapsize, centromere=gapsize, \
                contig=gapsize, clone=50000, \
                fragment=5000)

    def __init__(self, filename):

        super(Certificate, self).__init__(filename)

        fp = open(filename)
        self.lines = [CertificateLine(x) for x in fp.readlines()]

    def write(self, filename):
        fw = must_open(filename, "w")
        for b in self.lines:
            print >> fw, b

    def get_agp_gap(self, gap_type="contig"):
        gap_length = Certificate.gaps[gap_type]
        linkage = "yes" if gap_type in ("fragment", "clone") else "no"

        return ["N", gap_length, gap_type, linkage, ""]

    def write_AGP(self, filename, orientationguide={}, reindex=True):
        """
        For each component, we have two overlaps: North and South.

        =======
           ||||             South
           ====(=================)  Current BAC
           North             ||||
                             ===============

        For the case that says "Non-terminal", the overlap will not be
        considered. North-South would suggest a '+' orientation, South-North
        would suggest a '-' orientation. In most cases, unless the overlap
        involves phase1 BAC, the selected range will be shown as the brackets
        above - exclude North overlap, and include South overlap (aka the
        "left-greedy" rule).
        """
        fw = must_open(filename, "w")
        for aid, bb in groupby(self.lines, key=lambda x: x.aid):
            bb = list(bb)
            north, south = bb
            aid = north.aid
            assert aid == south.aid

            aphase = north.aphase
            chr = north.chr
            size = north.asize
            ar = [chr, 0, 0, 0]

            northline = southline = None
            northrange = southrange = None

            # Most gaps, except telomeres occur twice, so only do the "North"
            if north.is_gap:
                bar = ar + self.get_agp_gap(north.bid)
                northline = "\t".join(str(x) for x in bar)
            else:
                if north.isTerminal():
                    northrange = north.astart, north.astop

            if south.is_gap:
                if south.bid == "telomere":
                    bar = ar + self.get_agp_gap(south.bid)
                    southline = "\t".join(str(x) for x in bar)
            else:
                if south.isTerminal():
                    southrange = south.astart, south.astop
                else:
                    bar = ar + self.get_agp_gap("fragment")
                    southline = "\t".join(str(x) for x in bar)

            # Determine the orientation and clear range for the current BAC
            clr = [1, size]
            orientation = sorientation = None
            if northrange:
                start, stop = northrange
                Lhang = start - 1
                Rhang = size - stop

                orientation = '+' if Lhang < Rhang else '-'
                if north.bphase == 1 and north.bphase < aphase:
                    if Lhang < Rhang:  # North overlap at 5`
                        clr[0] = start
                    else:
                        clr[1] = stop
                # Override left-greedy (also see below)
                else:
                    if Lhang < Rhang:
                        clr[0] = stop + 1
                    else:
                        clr[1] = start - 1

            if southrange:
                start, stop = southrange
                Lhang = start - 1
                Rhang = size - stop

                sorientation = '+' if Lhang > Rhang else '-'
                # Override left-greedy (also see above)
                if aphase == 1 and aphase < south.bphase:
                    if Lhang < Rhang:  # South overlap at 5`
                        clr[0] = stop + 1
                    else:
                        clr[1] = start - 1
                else:
                    if Lhang < Rhang:
                        clr[0] = start
                    else:
                        clr[1] = stop

            if orientation:
                if sorientation:
                    try:
                        assert orientation == sorientation, \
                                "\n{0}\n{1}".format(north, south)
                    except AssertionError as e:
                        logging.debug(e)
            else:
                if sorientation:
                    orientation = sorientation
                else:  # Both overlaps fail to define orientation
                    orientation = orientationguide.get(aid, "+")

            component_type = "D" if aphase in (1, 2) else "F"
            bar = ar + [component_type, aid, clr[0], clr[1], orientation]
            cline = "\t".join(str(x) for x in bar)

            if northline:
                print >> fw, northline
            print >> fw, cline
            if southline:
                print >> fw, southline

        fw.close()

        if reindex:
            from jcvi.formats.agp import reindex
            reindex([filename])
            newagpfile = filename.replace(".agp", ".reindexed.agp")
            shutil.move(newagpfile, filename)


def main():

    actions = (
        ('bes', 'confirm the BES mapping'),
        ('flip', 'flip the FASTA sequences according to a set of references'),
        ('overlap', 'check terminal overlaps between two records'),
        ('neighbor', 'check neighbors of a component in agpfile'),
        ('blast', 'blast a component to componentpool'),
        ('certificate', 'make certificates for all overlaps in agpfile'),
        ('agp', 'make agpfile based on certificates'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def blast(args):
    """
    %prog blast allfasta clonename

    Insert a component into agpfile by aligning to the best hit in pool and see
    if they have good overlaps.
    """
    from jcvi.apps.command import run_megablast

    p = OptionParser(blast.__doc__)
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
    from jcvi.formats.blast import chain_HSPs

    p = OptionParser(overlap.__doc__)
    p.add_option("--dir", default=os.getcwd(),
            help="Download sequences to dir [default: %default]")
    p.add_option("--qreverse", default=False, action="store_true",
            help="Reverse seq a [default: %default]")
    p.add_option("--nochain", default=False, action="store_true",
            help="Do not chain adjacent HSPs [default: chain HSPs]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    afasta, bfasta = args
    dir = opts.dir
    chain = not opts.nochain

    # Check first whether it is file or accession name
    if not op.exists(afasta):
        afasta = op.join(dir, afasta + ".fasta")
        if not op.exists(afasta):  # Check to avoid redownload
            fetch([afasta, "--skipcheck", "--outdir=" + dir])

    if not op.exists(bfasta):
        bfasta = op.join(dir, bfasta + ".fasta")
        if not op.exists(bfasta):
            fetch([bfasta, "--skipcheck", "--outdir=" + dir])

    cmd = BLPATH("blastn")
    cmd += " -query {0} -subject {1}".format(afasta, bfasta)
    cmd += " -evalue 0.01 -outfmt 6 -perc_identity 99"

    fp = popen(cmd)
    hsps = fp.readlines()
    if len(hsps) == 0:
        print >> sys.stderr, "No match found."
        return None

    hsps = [BlastLine(x) for x in hsps]
    if chain:
        logging.debug("Chain HSPs in the Blast output.")
        hsps = chain_HSPs(hsps)

    besthsp = hsps[0]

    aid, asize = Fasta(afasta).itersizes().next()
    bid, bsize = Fasta(bfasta).itersizes().next()
    o = Overlap(besthsp, asize, bsize)
    o.print_graphic(qreverse=opts.qreverse)
    print >> sys.stderr, str(o)

    return o


@memoized
def phase(accession):
    gbdir = "gb"
    gbfile = op.join(gbdir, accession + ".gb")
    if not op.exists(gbfile):
        fetch([accession, "--skipcheck", "--outdir=" + gbdir, \
               "--format=gb"])
    rec = SeqIO.parse(gbfile, "gb").next()
    ph, keywords = get_phase(rec)
    return ph, len(rec)


def certificate(args):
    """
    %prog certificate tpffile certificatefile

    Generate certificate file for all overlaps in tpffile. tpffile can be
    generated by jcvi.formats.agp.tpf().

    North  chr1  2  0  AC229737.8  telomere     58443
    South  chr1  2  1  AC229737.8  AC202463.29  58443  37835  58443  + Non-terminal

    Each line describes a relationship between the current BAC and the
    north/south BAC. First, "North/South" tag, then the chromosome, phases of
    the two BACs, ids of the two BACs, the size and the overlap start-stop of
    the CURRENT BAC, and orientation. Each BAC will have two lines in the
    certificate file.
    """
    p = OptionParser(certificate.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    tpffile, certificatefile = args
    fastadir = "fasta"

    tpf = TPF(tpffile)

    data = {}
    if op.exists(certificatefile):
        # This will make updates resume-able and backed-up
        certificatefilebak = certificatefile + ".orig"
        shutil.copy2(certificatefile, certificatefilebak)

        fp = open(certificatefile)
        for row in fp:
            atoms = row.split()
            tag, aid, bid = atoms[0], atoms[4], atoms[5]
            data[(tag, aid, bid)] = row.strip()

    fw = must_open(certificatefile, "w")
    for i, a in enumerate(tpf):
        if a.is_gap:
            continue

        north, south = tpf.getNorthSouthClone(i)
        aid = a.component_id

        aphase, asize = phase(aid)

        for tag, p in (("North", north), ("South", south)):
            if not p:  # end of the chromosome
                ov = "telomere\t{0}".format(asize)
            elif p.isCloneGap:
                bphase = "0"
                ov = "{0}\t{1}".format(p.gap_type, asize)
            else:
                bid = p.component_id
                bphase, bsize = phase(bid)
                key = (tag, aid, bid)
                if key in data:
                    print >> fw,  data[key]
                    continue

                ar = [aid, bid, "--dir=" + fastadir]
                o = overlap(ar)
                ov = o.certificateline if o \
                        else "{0}\t{1}\tNone".format(bid, asize)

            print >> fw, "\t".join(str(x) for x in \
                    (tag, a.object, aphase, bphase, aid, ov))
            fw.flush()


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
    if not componentID in aorder:
        print >> sys.stderr, "Record {0} not present in `{1}`."\
                .format(componentID, agpfile)
        return

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


def agp(args):
    """
    %prog agp tpffile certificatefile agpfile

    Build agpfile from overlap certificates.

    Tiling Path File (tpf) is a file that lists the component and the gaps.
    It is a three-column file similar to belwo, also see jcvi.formats.agp.tpf():

    telomere        chr1
    AC229737.8      chr1    +
    AC202463.29     chr1    +

    Note: the orientation of the component is only used as a guide. If the
    orientation is derivable from a terminal overlap, it will use it regardless
    of what the tpf says.

    See jcvi.assembly.goldenpath.certificate() which generates a list of
    certificates based on agpfile. At first, it seems counter-productive to
    convert first agp to certificates then certificates back to agp.

    The certificates provide a way to edit the overlap information, so that the
    agpfile can be corrected (without changing agpfile directly).
    """
    from jcvi.formats.base import DictFile

    p = OptionParser(agp.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    tpffile, certificatefile, agpfile = args
    orientationguide = DictFile(tpffile, valuepos=2)
    cert = Certificate(certificatefile)
    cert.write_AGP(agpfile, orientationguide=orientationguide)


if __name__ == '__main__':
    main()
