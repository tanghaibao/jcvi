#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
From synteny blocks, reconstruct ancestral order by interleaving the genes in
between the anchors.
"""

import sys

from math import sqrt
from optparse import OptionParser

from jcvi.algorithms.synteny import AnchorFile, add_beds, check_beds
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('collinear', 'reduce synteny blocks to strictly collinear'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


# Non-linear transformation of anchor scores
score_convert = lambda x: int(sqrt(x))


def get_collinear(block):
    # block contains (gene a, gene b, score)
    asc_score, asc_chain = print_chain(block)
    desc_score, desc_chain = print_chain(block, ascending=False)
    return asc_chain if asc_score > desc_score else desc_chain


def print_chain(block, ascending=True):

    scope = 50  # reduce search complexity
    if not ascending:
        block = [(a, -b, c) for (a, b, c) in block]

    block.sort()
    bsize = len(block)
    fromm = [-1] * bsize
    scores = [score_convert(c) for (a, b, c) in block]

    for i, (a, b, c) in enumerate(block):
        for j in xrange(i + 1, i + scope):
            if j >= bsize:
                break

            d, e, f = block[j]

            # Ensure strictly collinear
            if d == a or b >= e:
                continue

            this_score = scores[i] + score_convert(f)
            if this_score > scores[j]:
                fromm[j] = i
                scores[j] = this_score

    scoresfromm = zip(scores, fromm)
    maxchain = max(scoresfromm)
    chainscore, chainend = maxchain
    solution = [scoresfromm.index(maxchain), chainend]
    last = chainend
    while True:
        _last = fromm[last]
        if _last == -1:
            break
        last = _last
        solution.append(last)

    solution.reverse()
    solution = [block[x] for x in solution]
    if not ascending:
        solution = [(a, -b, c) for (a, b, c) in solution]
    return chainscore, solution


def collinear(args):
    """
    %prog collinear a.b.anchors

    Reduce synteny blocks to strictly collinear, use dynamic programming in a
    procedure similar to DAGchainer.
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

        block = get_collinear(iblock)

        for q, s, score in block:
            q = qbed[q].accn
            s = sbed[s].accn
            print >> fw, "\t".join((q, s, str(score)))

    fw.close()


if __name__ == '__main__':
    main()
