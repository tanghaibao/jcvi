#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog rice.sorghum.last --qbed=rice.bed --sbed=sorghum.bed

Given a blast, we find the syntenic regions for every single gene. The
algorithm works by expanding the query gene to a window centered on the gene. A
single linkage algorithm follows that outputs the synteny block.

The result looks like the following:
Os01g0698300    Sb03g032090     S    7     +
Os01g0698500    Sb03g032140     G    11    +

The pairs (A, B) -- A is query, and then B is the syntenic region found.
G is "Gray gene", which means it does not have match to the region (fractionated
or inserted). In this case, a right flanker is used to represent the region.
S is "Syntelog", which means it has a match to the region. In this case, the match
itself is used to represent the region.  The number in the 4th column is the
synteny score. For the same query, it is ordered with decreasing synteny score.
The last column means orientation. "+" is same direction.
"""
from __future__ import print_function

import os.path as op
import logging
import sys
import sqlite3

from bisect import bisect_left
from itertools import groupby, tee

from jcvi.algorithms.lis import longest_increasing_subsequence, \
    longest_decreasing_subsequence
from jcvi.compara.synteny import check_beds, read_blast
from jcvi.utils.grouper import Grouper
from jcvi.formats.base import must_open
from jcvi.apps.base import OptionParser, OptionGroup


def transposed(data):
    x, y = zip(*data)
    return zip(y, x)


def get_flanker(group, query):
    """
    >>> get_flanker([(370, 15184), (372, 15178), (373, 15176), (400, 15193)],  385)
    ((373, 15176), (400, 15193), True)

    >>> get_flanker([(124, 13639), (137, 13625)], 138)
    ((137, 13625), (137, 13625), False)
    """
    group.sort()
    pos = bisect_left(group, (query, 0))
    left_flanker = group[0] if pos == 0 else group[pos-1]
    right_flanker = group[-1] if pos == len(group) else group[pos]
    # pick the closest flanker
    if abs(query - left_flanker[0]) < abs(query - right_flanker[0]):
        flanker, other = left_flanker, right_flanker
    else:
        flanker, other = right_flanker, left_flanker

    flanked = not (pos == 0 or pos == len(group) or flanker == query)

    return flanker, other, flanked


def find_synteny_region(query, sbed, data, window, cutoff, colinear=False):
    """
    Get all synteny blocks for a query, algorithm is single linkage
    anchors are a window centered on query

    Two categories of syntenic regions depending on what query is:
    (Syntelog): syntenic region is denoted by the syntelog
    (Gray gene): syntenic region is marked by the closest flanker
    """
    regions = []
    ysorted = sorted(data, key=lambda x: x[1])
    g = Grouper()

    a, b = tee(ysorted)
    next(b, None)
    for ia, ib in zip(a, b):
        pos1, pos2 = ia[1], ib[1]
        if pos2 - pos1 < window and sbed[pos1].seqid == sbed[pos2].seqid:
            g.join(ia, ib)

    for group in sorted(g):
        (qflanker, syntelog), (far_flanker, far_syntelog), flanked = \
            get_flanker(group, query)

        # run a mini-dagchainer here, take the direction that gives us most anchors
        if colinear:
            y_indexed_group = [(y, i) for i, (x, y) in enumerate(group)]
            lis = longest_increasing_subsequence(y_indexed_group)
            lds = longest_decreasing_subsequence(y_indexed_group)

            if len(lis) >= len(lds):
                track = lis
                orientation = "+"
            else:
                track = lds
                orientation = "-"

            group = [group[i] for (y, i) in track]

        xpos, ypos = zip(*group)
        score = min(len(set(xpos)), len(set(ypos)))

        if qflanker == query:
            gray = "S"
        else:
            gray = "G" if not flanked else "F"
            score -= 1  # slight penalty for not finding syntelog

        if score < cutoff:
            continue

        # y-boundary of the block
        left, right = group[0][1], group[-1][1]
        # this characterizes a syntenic region (left, right).
        # syntelog is -1 if it's a gray gene
        syn_region = (syntelog, far_syntelog, left,
                      right, gray, orientation, score)
        regions.append(syn_region)

    return sorted(regions, key=lambda x: -x[-1])  # decreasing synteny score


def batch_query(qbed, sbed, all_data, opts, fw=None, c=None, transpose=False):

    cutoff = int(opts.cutoff * opts.window)
    window = opts.window / 2
    colinear = opts.scoring == "collinear"
    qnote, snote = opts.qnote, opts.snote
    if qnote == "null" or snote == "null":
        qnote = op.basename(qbed.filename).split(".")[0]
        snote = op.basename(sbed.filename).split(".")[0]

    # process all genes present in the bed file
    if transpose:
        all_data = transposed(all_data)
        qbed, sbed = sbed, qbed
        qnote, snote = snote, qnote

    all_data.sort()
    def simple_bed(x): return (sbed[x].seqid, sbed[x].start)
    qsimplebed = qbed.simple_bed

    for seqid, ranks in groupby(qsimplebed, key=lambda x: x[0]):
        ranks = [x[1] for x in ranks]
        for r in ranks:
            rmin = max(r - window, ranks[0])
            rmax = min(r + window + 1, ranks[-1])
            rmin_pos = bisect_left(all_data, (rmin, 0))
            rmax_pos = bisect_left(all_data, (rmax, 0))
            data = all_data[rmin_pos:rmax_pos]
            regions = find_synteny_region(r, sbed, data, window,
                                          cutoff, colinear=colinear)
            for syntelog, far_syntelog, left, right, gray, orientation, score in regions:
                query = qbed[r].accn

                left_chr, left_pos = simple_bed(left)
                right_chr, right_pos = simple_bed(right)

                anchor = sbed[syntelog].accn
                anchor_chr, anchor_pos = simple_bed(syntelog)
                # below is useful for generating the syntenic region in the coge url
                left_dist = abs(anchor_pos - left_pos) \
                    if anchor_chr == left_chr else 0
                right_dist = abs(anchor_pos - right_pos) \
                    if anchor_chr == right_chr else 0
                flank_dist = (max(left_dist, right_dist) / 10000 + 1) * 10000

                far_syntelog = sbed[far_syntelog].accn

                left_pos, right_pos = sorted((left_pos, right_pos))
                data = [query, anchor, gray, score,
                        flank_dist, orientation, far_syntelog]
                pdata = data[:6] + [qnote, snote]
                if fw:
                    print("\t".join(str(x) for x in pdata), file=fw)
                    continue
                c.execute("insert into synteny values (?,?,?,?,?,?,?,?)", pdata)


def main(blastfile, p, opts):

    sqlite = opts.sqlite
    qbed, sbed, qorder, sorder, is_self = check_beds(blastfile, p, opts)
    filtered_blast = read_blast(blastfile, qorder, sorder,
                                is_self=is_self, ostrip=opts.strip_names)
    all_data = [(b.qi, b.si) for b in filtered_blast]

    c = None
    if sqlite:
        conn = sqlite3.connect(sqlite)
        c = conn.cursor()
        c.execute("drop table if exists synteny")
        c.execute("create table synteny (query text, anchor text, "
                  "gray varchar(1), score integer, dr integer, "
                  "orientation varchar(1), qnote text, snote text)")
        fw = None
    else:
        fw = must_open(opts.outfile, "w")

    batch_query(qbed, sbed, all_data, opts, fw=fw, c=c, transpose=False)
    if qbed.filename == sbed.filename:
        logging.debug("Self comparisons, mirror ignored")
    else:
        batch_query(qbed, sbed, all_data, opts, fw=fw, c=c, transpose=True)

    if sqlite:
        c.execute("create index q on synteny (query)")
        conn.commit()
        c.close()
    else:
        fw.close()


if __name__ == '__main__':

    p = OptionParser(__doc__)
    p.set_beds()
    p.set_stripnames()
    p.set_outfile()

    coge_group = OptionGroup(p, "CoGe-specific options")
    coge_group.add_option("--sqlite", help="Write sqlite database")
    coge_group.add_option("--qnote", default="null",
                          help="Query dataset group id")
    coge_group.add_option("--snote", default="null",
                          help="Subject dataset group id")

    params_group = OptionGroup(p, "Synteny parameters")
    params_group.add_option("--window", type="int", default=40,
                            help="Synteny window size")
    params_group.add_option("--cutoff", type="float", default=.1,
                            help="Minimum number of anchors to call synteny")
    supported_scoring = ("collinear", "density")
    params_group.add_option("--scoring", choices=supported_scoring,
                            default="collinear", help="Scoring scheme")

    p.add_option_group(coge_group)
    p.add_option_group(params_group)

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(not p.print_help())

    blastfile, = args
    main(blastfile, p, opts)
