#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Based on read pair mappings, construct contig graph
"""

import sys

from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.bed import BedLine
from jcvi.utils.iter import pairwise
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('bed', 'construct links from bed file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed bedfile

    Construct contig links based on bed file
    """
    p = OptionParser(bed.__doc__)
    p.add_option("-l", dest="links", type="int", default=3,
            help="Minimum number of mate pairs to bundle [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    bedfile, = args
    linksfile = bedfile.rsplit(".", 1)[0] + ".links"

    fp = open(bedfile)
    
    contigGraph = defaultdict(int)
    for a, b in pairwise(fp):
        a, b = BedLine(a), BedLine(b)
        """
        Criteria for valid contig edge
        1. for/rev not mapping to the same scaffold
        2. TODO: check boundary and insert size
        """
        if a.accn[:-1] != b.accn[:-1]: continue
        aseqid, bseqid = a.seqid, b.seqid
        if aseqid == bseqid: continue

        adir = a.stuff[1] 
        bdir = b.stuff[1]

        if aseqid > bseqid:
            aseqid, bseqid = bseqid, aseqid

        contigGraph[(aseqid, bseqid)] += 1 

    for pair, mates in sorted(contigGraph.items()):
        if mates < opts.links: continue 
        print pair, mates


if __name__ == '__main__':
    main()
