#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Catalog gene losses, and bites within genes.
"""

import sys

from optparse import OptionParser
from itertools import groupby

from jcvi.formats.blast import Blast
from jcvi.utils.range import range_minmax, range_overlap
from jcvi.utils.cbook import gene_name
from jcvi.algorithms.synteny import add_beds, check_beds
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('loss', 'extract likely gene loss candidates'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def region_str(region):
    return "{0}:{1}-{2}".format(*region)


def loss(args):
    """
    %prog loss a.b.i1.blocks a.b-genomic.blast

    Extract likely gene loss candidates between genome a and b.
    """
    p = OptionParser(loss.__doc__)
    p.add_option("--gdist", default=20,
                 help="Gene distance [default: %default]")
    p.add_option("--bdist", default=20000,
                 help="Base pair distance [default: %default]")
    add_beds(p)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blocksfile, genomicblast = args
    gdist, bdist = opts.gdist, opts.bdist
    qbed, sbed, qorder, sorder, is_self = check_beds(blocksfile, p, opts)
    blocks = []
    fp = open(blocksfile)
    genetrack = {}
    proxytrack = {}
    for row in fp:
        a, b = row.split()
        genetrack[a] = b
        blocks.append((a, b))

    data = []
    for key, rows in groupby(blocks, key=lambda x: x[-1]):
        rows = list(rows)
        data.append((key, rows))

    imax = len(data) - 1
    for i, (key, rows) in enumerate(data):
        if i == 0 or i == imax:
            continue
        if key != '.':
            continue

        before, br = data[i - 1]
        after, ar = data[i + 1]
        bi, bx = sorder[before]
        ai, ax = sorder[after]
        dist = abs(bi - ai)
        if bx.seqid != ax.seqid or dist > gdist:
            continue

        start, end = range_minmax(((bx.start, bx.end), (ax.start, ax.end)))
        proxy = (bx.seqid, start - bdist, end + bdist)
        for a, b in rows:
            proxytrack[a] = proxy

    blast = Blast(genomicblast)
    tags = {}
    for query, bb in blast.iter_hits():
        query = gene_name(query)
        if query not in proxytrack:
            continue

        proxy = proxytrack[query]
        tag = "NS"
        for b in bb:
            hsp = (b.subject, b.sstart, b.sstop)
            if range_overlap(proxy, hsp):
                proxytrack[query] = hsp
                tag = "S"
                break
        tags[query] = tag

    for b in qbed:
        accn = b.accn
        target_region = genetrack[accn]
        if accn in proxytrack:
            target_region = region_str(proxytrack[accn])
            if accn in tags:
                target_region += "[{0}]".format(tags[accn])
            else:
                target_region += "[NF]"
        print "\t".join((accn, target_region))


if __name__ == '__main__':
    main()
