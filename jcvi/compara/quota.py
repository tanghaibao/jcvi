#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Quota synteny alignment (QUOTA-ALIGN)

%prog [options] anchorsfile --qbed=qbedfile --sbed=sbedfile

This python program does the following:
1. merge 2D-overlapping blocks (now skipped, but existed in original version)
2. build constraints that represent 1D-overlap among blocks
3. feed the data into the linear programming solver

The algorithm is described in Tang et al. BMC Bioinformatics 2011.
"Screening synteny blocks in pairwise genome comparisons through integer
programming."
"""
from __future__ import print_function

import os.path as op
import sys
from six.moves import StringIO
import logging

from jcvi.utils.range import range_overlap
from jcvi.utils.grouper import Grouper
from jcvi.algorithms.lpsolve import GLPKSolver, SCIPSolver
from jcvi.compara.synteny import AnchorFile, _score, check_beds
from jcvi.formats.base import must_open
from jcvi.apps.base import OptionParser


def get_1D_overlap(eclusters, depth=1):
    """
    Find blocks that are 1D overlapping,
    returns cliques of block ids that are in conflict
    """
    overlap_set = set()
    active = set()

    ends = []
    for i, (chr, left, right) in enumerate(eclusters):
        ends.append((chr, left, 0, i))  # 0/1 for left/right-ness
        ends.append((chr, right, 1, i))
    ends.sort()

    chr_last = ""
    for chr, pos, left_right, i in ends:
        if chr != chr_last:
            active.clear()
        if left_right == 0:
            active.add(i)
        else:
            active.remove(i)

        if len(active) > depth:
            overlap_set.add(tuple(sorted(active)))

        chr_last = chr

    return overlap_set


def get_2D_overlap(chain, eclusters):
    """
    Implements a sweep line algorithm, that has better running time than naive O(n^2):
    assume block has x_ends, and y_ends for the bounds

    1. sort x_ends, and take a sweep line to scan the x_ends
    2. if left end, test y-axis intersection of current block with `active` set;
       also put this block in the `active` set
    3. if right end, remove block from the `active` set
    """
    mergeables = Grouper()
    active = set()

    x_ends = []
    for i, (range_x, range_y, score) in enumerate(eclusters):
        chr, left, right = range_x
        x_ends.append((chr, left, 0, i))  # 0/1 for left/right-ness
        x_ends.append((chr, right, 1, i))
    x_ends.sort()

    chr_last = ""
    for chr, pos, left_right, i in x_ends:
        if chr != chr_last:
            active.clear()
        if left_right == 0:
            active.add(i)
            for x in active:
                # check y-overlap
                if range_overlap(eclusters[x][1], eclusters[i][1]):
                    mergeables.join(x, i)
        else:  # right end
            active.remove(i)

        chr_last = chr

    return mergeables


def make_range(clusters, extend=0):
    """
    Convert to interval ends from a list of anchors
    extend modifies the xmax, ymax boundary of the box,
    which can be positive or negative
    very useful when we want to make the range as fuzzy as we specify
    """
    eclusters = []
    for cluster in clusters:
        xlist, ylist, scores = zip(*cluster)
        score = _score(cluster)

        xchr, xmin = min(xlist)
        xchr, xmax = max(xlist)
        ychr, ymin = min(ylist)
        ychr, ymax = max(ylist)

        # allow fuzziness to the boundary
        xmax += extend
        ymax += extend
        # because extend can be negative values, we don't want it to be less than min
        if xmax < xmin:
            xmin, xmax = xmax, xmin
        if ymax < ymin:
            ymin, ymax = ymax, ymin

        eclusters.append(((xchr, xmin, xmax), (ychr, ymin, ymax), score))

    return eclusters


def get_constraints(clusters, quota=(1, 1), Nmax=0):
    """
    Check pairwise cluster comparison, if they overlap then mark edge as conflict
    """
    qa, qb = quota
    eclusters = make_range(clusters, extend=-Nmax)
    # (1-based index, cluster score)
    nodes = [(i + 1, c[-1]) for i, c in enumerate(eclusters)]

    eclusters_x, eclusters_y, scores = zip(*eclusters)

    # represents the contraints over x-axis and y-axis
    constraints_x = get_1D_overlap(eclusters_x, qa)
    constraints_y = get_1D_overlap(eclusters_y, qb)

    return nodes, constraints_x, constraints_y


def format_lp(nodes, constraints_x, qa, constraints_y, qb):
    """
    Maximize
     4 x1 + 2 x2 + 3 x3 + x4
    Subject To
     x1 + x2 <= 1
    End
    """
    lp_handle = StringIO()

    lp_handle.write("Maximize\n ")
    records = 0
    for i, score in nodes:
        lp_handle.write("+ %d x%d " % (score, i))
        # SCIP does not like really long string per row
        records += 1
        if records % 10 == 0:
            lp_handle.write("\n")
    lp_handle.write("\n")

    num_of_constraints = 0
    lp_handle.write("Subject To\n")
    for c in constraints_x:
        additions = " + ".join("x%d" % (x + 1) for x in c)
        lp_handle.write(" %s <= %d\n" % (additions, qa))
    num_of_constraints += len(constraints_x)

    # non-self
    if not (constraints_x is constraints_y):
        for c in constraints_y:
            additions = " + ".join("x%d" % (x + 1) for x in c)
            lp_handle.write(" %s <= %d\n" % (additions, qb))
        num_of_constraints += len(constraints_y)

    print(
        "number of variables (%d), number of constraints (%d)"
        % (len(nodes), num_of_constraints),
        file=sys.stderr,
    )

    lp_handle.write("Binary\n")
    for i, score in nodes:
        lp_handle.write(" x%d\n" % i)

    lp_handle.write("End\n")

    lp_data = lp_handle.getvalue()
    lp_handle.close()

    return lp_data


def solve_lp(
    clusters,
    quota,
    work_dir="work",
    Nmax=0,
    self_match=False,
    solver="SCIP",
    verbose=False,
):
    """
    Solve the formatted LP instance
    """
    qb, qa = quota  # flip it
    nodes, constraints_x, constraints_y = get_constraints(clusters, (qa, qb), Nmax=Nmax)

    if self_match:
        constraints_x = constraints_y = constraints_x | constraints_y

    lp_data = format_lp(nodes, constraints_x, qa, constraints_y, qb)

    if solver == "SCIP":
        filtered_list = SCIPSolver(lp_data, work_dir, verbose=verbose).results
        if not filtered_list:
            print("SCIP fails... trying GLPK", file=sys.stderr)
            filtered_list = GLPKSolver(lp_data, work_dir, verbose=verbose).results

    elif solver == "GLPK":
        filtered_list = GLPKSolver(lp_data, work_dir, verbose=verbose).results
        if not filtered_list:
            print("GLPK fails... trying SCIP", file=sys.stderr)
            filtered_list = SCIPSolver(lp_data, work_dir, verbose=verbose).results

    return filtered_list


def read_clusters(qa_file, qorder, sorder):
    af = AnchorFile(qa_file)
    blocks = af.blocks
    clusters = []
    for block in blocks:
        cluster = []
        for a, b, score in block:
            ia, oa = qorder[a]
            ib, ob = sorder[b]
            ca, cb = oa.seqid, ob.seqid
            cluster.append(((ca, ia), (cb, ib), score))
        clusters.append(cluster)

    return clusters


def main(args):
    p = OptionParser(__doc__)

    p.set_beds()
    p.add_option(
        "--quota",
        default="1:1",
        help="`quota mapping` procedure -- screen blocks to constrain mapping"
        " (useful for orthology), "
        "put in the format like (#subgenomes expected for genome X):"
        "(#subgenomes expected for genome Y) "
        "[default: %default]",
    )
    p.add_option(
        "--Nm",
        dest="Nmax",
        type="int",
        default=10,
        help="distance cutoff to tolerate two blocks that are "
        "slightly overlapping (cutoff for `quota mapping`) "
        "[default: %default units (gene or bp dist)]",
    )

    supported_solvers = ("SCIP", "GLPK")
    p.add_option(
        "--self",
        dest="self_match",
        action="store_true",
        default=False,
        help="you might turn this on when screening paralogous blocks, "
        "esp. if you have reduced mirrored blocks into non-redundant set",
    )
    p.add_option(
        "--solver",
        default="SCIP",
        choices=supported_solvers,
        help="use MIP solver [default: %default]",
    )
    p.set_verbose(help="Show verbose solver output")

    p.add_option(
        "--screen",
        default=False,
        action="store_true",
        help="generate new anchors file [default: %default]",
    )

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    qa_file, = args
    qbed, sbed, qorder, sorder, is_self = check_beds(qa_file, p, opts)

    # sanity check for the quota
    if opts.quota:
        try:
            qa, qb = opts.quota.split(":")
            qa, qb = int(qa), int(qb)
        except:
            print(
                "quota string should be the form x:x (2:4, 1:3, etc.)", file=sys.stderr
            )
            sys.exit(1)

        if opts.self_match and qa != qb:
            raise Exception(
                "when comparing genome to itself, "
                "quota must be the same number "
                "(like 1:1, 2:2) you have %s" % opts.quota
            )
        quota = (qa, qb)

    self_match = opts.self_match

    clusters = read_clusters(qa_file, qorder, sorder)
    for cluster in clusters:
        assert len(cluster) > 0

    # below runs `quota mapping`
    work_dir = op.join(op.dirname(op.abspath(qa_file)), "work")

    selected_ids = solve_lp(
        clusters,
        quota,
        work_dir=work_dir,
        Nmax=opts.Nmax,
        self_match=self_match,
        solver=opts.solver,
        verbose=opts.verbose,
    )

    logging.debug("Selected {0} blocks.".format(len(selected_ids)))
    prefix = qa_file.rsplit(".", 1)[0]
    suffix = "{0}x{1}".format(qa, qb)
    outfile = ".".join((prefix, suffix))
    fw = must_open(outfile, "w")
    print(",".join(str(x) for x in selected_ids), file=fw)
    fw.close()
    logging.debug("Screened blocks ids written to `{0}`.".format(outfile))

    if opts.screen:
        from jcvi.compara.synteny import screen

        new_qa_file = ".".join((prefix, suffix, "anchors"))
        largs = [qa_file, new_qa_file, "--ids", outfile]
        if opts.qbed and opts.sbed:
            largs += ["--qbed={0}".format(opts.qbed)]
            largs += ["--sbed={0}".format(opts.sbed)]
        screen(largs)


if __name__ == "__main__":
    main(sys.argv[1:])
