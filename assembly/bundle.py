#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Based on read pair mappings, construct contig graph
"""

import sys
import logging

from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.bed import BedLine
from jcvi.formats.sizes import Sizes
from jcvi.formats.blast import report_pairs
from jcvi.utils.iter import pairwise
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('bed', 'construct links from bed file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def calculate_hang(match_start, match_end, strand, size):
    if strand=='+':
        return size - match_start
    else:
        return match_end


def bed(args):
    """
    prog bed bedfile fastafile

    Construct contig links based on bed file
    """
    p = OptionParser(bed.__doc__)
    p.add_option("-l", dest="links", type="int", default=2,
            help="Minimum number of mate pairs to bundle [default: %default]")
    p.add_option("--debug", dest="debug", default=False, action="store_true",
            help="Print verbose debug infor about checking mates [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    bedfile, fastafile = args
    links = opts.links
    debug = opts.debug

    linksfile = bedfile.rsplit(".", 1)[0] + ".links"

    sizes = Sizes(fastafile)

    fp = open(bedfile)
    
    sample = 100000
    logging.debug(\
        "Sample only the first {0} lines to estimate inserts".format(sample))
    beds = [BedLine(x) for i, x in enumerate(fp) if i < sample]
    meandist, stdev = report_pairs(beds, cutoff=0, dialect="bed")
    hangs_cutoff = int(meandist + 1.96 * stdev)
    logging.debug("Size of hangs have to be <= {0} bases.".format(hangs_cutoff))

    contigGraph = defaultdict(list)
    rs = lambda x: x.accn[:-1]

    fp.seek(0)

    for a, b in pairwise(fp):
        """
        Criteria for valid contig edge
        1. for/rev do not mapping to the same scaffold (useful for linking)
        2. assuming innie (outie has to be flipped first), order the contig pair
        3. calculate sequence hangs, valid hangs are smaller than insert size
        """
        a, b = BedLine(a), BedLine(b)
        if debug:
            print >> sys.stderr, a
            print >> sys.stderr, b

        if rs(a) != rs(b): continue
        pe = rs(a)

        if a.seqid == b.seqid: continue
        """
        The algorithm that determines the oNo of this contig pair, the contig
        order is determined by +-, assuming that the matepairs are `innies`. In
        below, we determine the order and orientation by flipping if necessary,
        bringing the strandness of two contigs to the expected +-.
        """
        if a.strand == b.strand:
            if b.strand == "+": # case: "++", flip b
                b.reverse_complement(sizes)
            else: # case: "--", flip a
                a.reverse_complement(sizes)

        if b.strand == "+": # case: "-+"
            a, b = b, a

        assert a.strand == "+" and b.strand == "-"
        """
        ------===----          -----====----
              |_ahang            bhang_| 
        """
        aseqid = a.seqid[:-1] if a.seqid[-1]=='-' else a.seqid
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
        if hangs > meandist + 1.96 * stdev: 
            if debug:
                print >> sys.stderr, "invalid link. skipped." 
            continue

        # Ignore redundant mates, in this case if the size of the hangs is seen,
        # then it is probably redundant, since indepenet circularization event
        # will give slightly different value)
        if hangs in [x[1] for x in contigGraph[(a.seqid, b.seqid)]]: continue
        contigGraph[(a.seqid, b.seqid)].append((pe, hangs))

    for pair, mates in sorted(contigGraph.items()):
        if len(mates) < links: continue 
        print pair, mates


if __name__ == '__main__':
    main()
