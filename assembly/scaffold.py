#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Yet another tool to perform scaffolding based on long reads or alternative assembly.

Note several other modules already exist: sopra, bambus, syntenypath that each
does scaffolding differently.
"""

import sys

from itertools import groupby

from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('connect', 'connect contigs using long reads'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def connect(args):
    """
    %prog connect assembly.fasta read_mapping.blast

    Connect contigs using long reads.
    """
    from jcvi.formats.sizes import Sizes
    from jcvi.formats.blast import Blast
    from jcvi.utils.iter import pairwise
    from jcvi.utils.range import range_intersect
    from jcvi.algorithms.graph import BiGraph, BiEdge
    from jcvi.assembly.syntenypath import graph_to_agp

    p = OptionParser(connect.__doc__)
    p.add_option("--clip", default=2000, type="int",
            help="Only consider end of contigs [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, blastfile = args
    clip = opts.clip

    sizes = Sizes(fastafile).mapping
    blast = Blast(blastfile)
    blasts = []
    for b in blast:
        seqid = b.subject
        size = sizes[seqid]
        start, end = b.sstart, b.sstop
        cstart, cend = min(size, clip), max(0, size - clip)
        if start > cstart and end < cend:
            continue
        blasts.append(b)

    key = lambda x: x.query
    blasts.sort(key=key)
    g = BiGraph()
    for query, bb in groupby(blasts, key=key):
        bb = sorted(bb, key=lambda x: x.qstart)
        nsubjects = len(set(x.subject for x in bb))
        if nsubjects == 1:
            continue
        print "\n".join(str(x) for x in bb)
        for a, b in pairwise(bb):
            astart, astop = a.qstart, a.qstop
            bstart, bstop = b.qstart, b.qstop
            if a.subject == b.subject:
                continue

            arange = astart, astop
            brange = bstart, bstop
            ov = range_intersect(arange, brange)
            alen = astop - astart + 1
            blen = bstop - bstart + 1
            if ov:
                ostart, ostop = ov
                ov = ostop - ostart + 1

            print ov, alen, blen
            if ov and (ov > alen / 2 or ov > blen / 2):
                print "Too much overlap ({0})".format(ov)
                continue

            asub = a.subject
            bsub = b.subject
            atag = ">" if a.orientation == "+" else "<"
            btag = ">" if b.orientation == "+" else "<"
            e = BiEdge(asub, bsub, atag, btag)
            g.add_edge(e)
            print "=" * 5, e

    graph_to_agp(g, blastfile, fastafile, verbose=False)


if __name__ == '__main__':
    main()
