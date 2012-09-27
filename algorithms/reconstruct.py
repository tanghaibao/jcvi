#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
From synteny blocks, reconstruct ancestral order by interleaving the genes in
between the anchors.
"""

import sys

from optparse import OptionParser

from jcvi.algorithms.synteny import AnchorFile, add_beds, check_beds
from jcvi.algorithms.lis import heaviest_increasing_subsequence
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('collinear', 'reduce synteny blocks to strictly collinear'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_collinear(block):
    # block contains (gene a, gene b, score)
    iblock = list(enumerate(block))
    iblock_subsorted = sorted(iblock, key=lambda (i, x): x)
    ib = [(i, x[-1]) for i, x in iblock_subsorted]
    his = heaviest_increasing_subsequence(ib)

    ii, bb = zip(*his)
    return ii


def collinear(args):
    """
    %prog collinear a.b.anchors

    Reduce synteny blocks to strictly collinear.
    # TODO: THIS IS NOT WORKING YET
    """
    p = OptionParser(collinear.__doc__)
    add_beds(p)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    anchorfile, = args
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)

    af = AnchorFile(anchorfile)
    newanchorfile = anchorfile.rsplit(".", 1)[0] + ".collinear.anchors"
    fw = open(newanchorfile, "w")

    blocks = af.blocks
    for block in blocks:
        print >> fw, "#" * 3
        iblock = []
        for q, s, score in block:
            qi, q = qorder[q]
            si, s = sorder[s]
            score = int(long(score))
            iblock.append([qi, si, score])

        ii = get_collinear(iblock)

        block = [block[x] for x in ii]
        for q, s, score in block:
            print >> fw, "\t".join((q, s, str(score)))

    fw.close()


if __name__ == '__main__':
    main()
