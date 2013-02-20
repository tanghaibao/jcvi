#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Catalog gene losses, and bites within genes.
"""

import sys

from optparse import OptionParser
from itertools import groupby

from jcvi.formats.blast import Blast
from jcvi.formats.bed import Bed
from jcvi.utils.range import range_minmax, range_overlap
from jcvi.utils.cbook import gene_name
from jcvi.algorithms.synteny import add_beds, check_beds
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('loss', 'extract likely gene loss candidates'),
        # Specific study (requires specific datasets)
        ('napus', 'extract napus gene loss vs diploid ancestors'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_tag(name, order):
    if name[0] == '[':
        tag, tname = name[1:].split(']')
        seqid, se = tname.split(":")
        start, end = se.split("-")
        start, end = int(start), int(end)
    else:
        tag = None
        xi, x = order[name]
        seqid, start, end = x.seqid, x.start, x.end
    return tag, (seqid, start, end)


def napus(args):
    """
    %prog napus napus.bed brapa.boleracea.i1.blocks diploid.napus.fractionation

    Extract napus gene loss vs diploid ancestors. We are looking specifically
    for anything that has the pattern:

        BR - BO    or     BR - BO
        |                       |
        AN                     CN

    Step 1: extract BR - BO syntenic pairs
    Step 2: get diploid gene retention patterns from BR or BO as query
    Step 3: look for if AN or CN is NS(non-syntenic) or NF(not found) and
    specifically with NS, the NS location is actually the homeologous site.
    """
    p = OptionParser(napus.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    napusbed, brbo, dpnp = args
    retention = {}
    fp = open(dpnp)
    for row in fp:
        seqid, query, hit = row.split()
        retention[query] = hit

    order = Bed(napusbed).order

    fp = open(brbo)
    for row in fp:
        br, bo = row.split()
        if '.' in (br, bo):
            continue
        an, cn = retention[br], retention[bo]
        if '.' in (an, cn):
            continue

        row = "\t".join((br, bo, an, cn))

        # label loss candidates
        antag, anrange = get_tag(an, order)
        cntag, cnrange = get_tag(cn, order)

        if range_overlap(anrange, cnrange):
            if (antag, cntag) == ("NS", None):
                row = row + "\tAN LOST"
            if (antag, cntag) == (None, "NS"):
                row = row + "\tCN LOST"

        print row


def region_str(region):
    return "{0}:{1}-{2}".format(*region)


def loss(args):
    """
    %prog loss a.b.i1.blocks a.b-genomic.blast

    Extract likely gene loss candidates between genome a and b.
    """
    p = OptionParser(loss.__doc__)
    p.add_option("--bed", default=False, action="store_true",
                 help="Genomic BLAST is in bed format [default: %default]")
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

    tags = {}
    if opts.bed:
        bed = Bed(genomicblast, sorted=False)
        key = lambda x: gene_name(x.accn.rsplit(".", 1)[0])
        for query, bb in groupby(bed, key=key):
            bb = list(bb)
            if query not in proxytrack:
                continue

            proxy = proxytrack[query]
            tag = "NS"
            best_b = bb[0]
            for b in bb:
                hsp = (b.seqid, b.start, b.end)
                if range_overlap(proxy, hsp):
                    tag = "S"
                    best_b = b
                    break

            hsp = (best_b.seqid, best_b.start, best_b.end)
            proxytrack[query] = hsp
            tags[query] = tag

    else:
        blast = Blast(genomicblast)
        for query, bb in blast.iter_hits():
            bb = list(bb)
            query = gene_name(query)
            if query not in proxytrack:
                continue

            proxy = proxytrack[query]
            tag = "NS"
            best_b = bb[0]
            for b in bb:
                hsp = (b.subject, b.sstart, b.sstop)
                if range_overlap(proxy, hsp):
                    tag = "S"
                    best_b = b
                    break

            hsp = (best_b.subject, best_b.sstart, best_b.sstop)
            proxytrack[query] = hsp
            tags[query] = tag

    for b in qbed:
        accn = b.accn
        target_region = genetrack[accn]
        if accn in proxytrack:
            target_region = region_str(proxytrack[accn])
            if accn in tags:
                ptag = "[{0}]".format(tags[accn])
            else:
                ptag = "[NF]"
            target_region = ptag + target_region

        print "\t".join((b.seqid, accn, target_region))


if __name__ == '__main__':
    main()
