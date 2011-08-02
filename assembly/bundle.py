#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Based on read pair mappings, construct contig graph
"""

import sys
import logging

from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.base import must_open
from jcvi.formats.bed import BedLine, pairs
from jcvi.formats.sizes import Sizes
from jcvi.utils.iter import pairwise
from jcvi.apps.base import ActionDispatcher, debug
debug()


class LinkLine (object):

    def __init__(self, row):
        args = row.split()
        # MKDU973T  mte1-26c10_002  mte1-26c10_011  +-  5066  23563
        self.mate = args[0]
        self.aseqid = args[1]
        self.bseqid = args[2]
        self.orientation = args[3]
        self.insert = int(args[4])
        self.distance = int(args[5])


class ContigLink (object):

    def __init__(self, a, b, insert=3000, cutoff=6000):
        assert isinstance(a, BedLine) and isinstance(b, BedLine)
        self.a = a
        self.b = b
        self.insert = insert
        assert self.insert > 0

        self.cutoff = cutoff
        assert self.cutoff > self.insert

    def __str__(self):
        aseqid = self.a.seqid.rstrip("-")
        bseqid = self.b.seqid.rstrip("-")
        return "\t".join(str(x) for x in (aseqid, bseqid,
            self.orientation, self.insert, self.distance))

    def get_orientation(self, aseqid, bseqid):
        """
        Determine N/A/I/O from the names.
        """
        ta = '-' if aseqid[-1] == '-' else '+'
        tb = '-' if bseqid[-1] == '-' else '+'
        pair = ta + tb
        assert pair in ('++', '--', '+-', '-+')

        return pair

    def flip_innie(self, sizes, debug=False):
        """
        The algorithm that determines the oNo of this contig pair, the contig
        order is determined by +-, assuming that the matepairs are `innies`. In
        below, we determine the order and orientation by flipping if necessary,
        bringing the strandness of two contigs to the expected +-.

        sizes: the contig length dictionary Sizes().
        """
        a, b = self.a, self.b
        Pa, Pb = a.start - 1, b.start - 1
        Ea, Eb = a.end, b.end

        if a.strand == b.strand:
            if b.strand == "+":  # case: "++", flip b
                b.reverse_complement(sizes)
            else:  # case: "--", flip a
                a.reverse_complement(sizes)

        if b.strand == "+":  # case: "-+"
            a, b = b, a
            Pa, Pb = Pb, Pa
            Ea, Eb = Eb, Ea

        assert a.strand == "+" and b.strand == "-"
        """
        ------===----          -----====----
              |_ahang            bhang_|
        """
        aseqid = a.seqid.rstrip('-')
        size = sizes.get_size(aseqid)
        ahang = size - a.start + 1
        bhang = b.end

        if debug:
            print >> sys.stderr, '*' * 60
            print >> sys.stderr, a
            print >> sys.stderr, b
            print >> sys.stderr, "ahang={0} bhang={1}".format(ahang, bhang)

        # Valid links need to be separated by the lib insert
        hangs = ahang + bhang
        if hangs > self.cutoff:
            if debug:
                print >> sys.stderr, "invalid link ({0}).".format(hangs)
            return False

        pair = self.get_orientation(a.seqid, b.seqid)

        insert = self.insert
        if pair == "++":    # Normal
            distance = insert + Pa - Eb
        elif pair == "--":  # Anti-normal
            distance = insert - Ea + Pb
        elif pair == "+-":  # Innie
            distance = insert + Pa + Pb
        elif pair == "-+":  # Outie
            distance = insert - Ea - Eb

        # Pair (1+, 2-) is the same as (2+, 1-), only do the canonical one
        if a.seqid > b.seqid:
            a.reverse_complement(sizes)
            b.reverse_complement(sizes)
            a, b = b, a

        self.a, self.b = a, b
        self.distance = distance
        self.orientation = self.get_orientation(a.seqid, b.seqid)

        return True


def main():

    actions = (
        ('link', 'construct links from bed file'),
        ('bundle', 'bundle multiple links into contig edges'),
        ('query', 'query the path between given contigs from the bundle'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def query(args):
    """
    %prog query bundlefile sourcectg targetctg

    Query the path from sourcectg to targetctg using links in the bundlefile.
    """
    from jcvi.algorithms.graph import nx, shortest_path

    p = OptionParser(query.__doc__)

    opts, args = p.parse_args(args)

    if len(args) < 3:
        sys.exit(not p.print_help())

    bundlefile, srcctg, targetctg = args

    g = nx.MultiGraph()
    fp = open(bundlefile)
    for row in fp:
        # 5  contig_12167  contig_43478  -+  0  -664
        c = LinkLine(row)
        g.add_edge(c.aseqid, c.bseqid)

    ctgs = shortest_path(g, srcctg, targetctg)

    ctgs = set(ctgs)
    fp = open(bundlefile)
    for row in fp:
        c = LinkLine(row)
        if c.aseqid in ctgs and c.bseqid in ctgs:
            print row.rstrip()


def bundle(args):
    """
    %prog bundle linkfiles

    Bundle contig links into high confidence contig edges. This is useful to
    combine multiple linkfiles (from different libraries).
    """
    import numpy as np

    p = OptionParser(bundle.__doc__)
    p.add_option("--links", type="int", default=2,
            help="Minimum number of mate pairs to bundle [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fp = must_open(args)
    contigGraph = defaultdict(list)
    for row in fp:
        c = LinkLine(row)
        contigGraph[(c.aseqid, c.bseqid, c.orientation)].append(c.distance)

    for (aseqid, bseqid, orientation), distances in contigGraph.items():
        mates = len(distances)
        if mates < opts.links:
            continue

        distance = int(np.median(distances))
        print "\t".join(str(x) for x in \
                (mates, aseqid, bseqid, orientation, 0, distance))


def link(args):
    """
    %prog link bedfile fastafile

    Construct contig links based on bed file. Use --prefix to limit the links
    between contigs that start with the same prefix_xxx.
    """
    p = OptionParser(link.__doc__)
    p.add_option("--insert", type="int", default=0,
            help="Mean insert size [default: estimate from data]")
    p.add_option("--cutoff", type="int", default=0,
            help="Largest distance expected for linkage " + \
                 "[default: estimate from data]")
    p.add_option("--prefix", default=False, action="store_true",
            help="Only keep links between IDs with same prefix [default: %default]")
    p.add_option("--debug", dest="debug", default=False, action="store_true",
            help="Print verbose info when checking mates [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, fastafile = args
    debug = opts.debug
    cutoff = opts.cutoff

    sizes = Sizes(fastafile)

    cutoffopt = "--cutoff={0}".format(cutoff)
    mateorientationopt = '--mateorientation=+-'
    bedfile, (meandist, stdev, p0, p1, p2) = \
            pairs([bedfile, cutoffopt, mateorientationopt])

    maxcutoff = cutoff or p2
    insert = opts.insert or p0
    logging.debug("Mate hangs must be <= {0}, --cutoff to override".\
            format(maxcutoff))

    rs = lambda x: x.accn[:-1]

    fp = open(bedfile)
    linksfile = bedfile.rsplit(".", 1)[0] + ".links"
    fw = open(linksfile, "w")

    for a, b in pairwise(fp):
        """
        Criteria for valid contig edge
        1. for/rev do not mapping to the same scaffold (useful for linking)
        2. assuming innie (outie must be flipped first), order the contig pair
        3. calculate sequence hangs, valid hangs are smaller than insert size
        """
        a, b = BedLine(a), BedLine(b)

        if rs(a) != rs(b):
            continue
        pe = rs(a)

        # Intra-contig links
        if a.seqid == b.seqid:
            continue

        # Use --prefix to limit the links between seqids with the same prefix
        # For example, contigs of the same BAC, mth2-23j10_001, mth-23j10_002
        if opts.prefix:
            aprefix = a.seqid.split("_")[0]
            bprefix = b.seqid.split("_")[0]
            if aprefix != bprefix:
                continue

        cl = ContigLink(a, b, insert=insert, cutoff=maxcutoff)
        if cl.flip_innie(sizes, debug=debug):
            print >> fw, "\t".join((pe, str(cl)))


if __name__ == '__main__':
    main()
