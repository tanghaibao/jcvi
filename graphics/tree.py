#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog newicktree

Plot Newick formatted tree.
"""


import sys
import logging

from optparse import OptionParser

from ete2 import Tree
from jcvi.graphics.base import plt, _, set_image_options
from jcvi.apps.base import debug
debug()


def draw_tree(ax, tx, rmargin=.2, outgroup=None):

    t = Tree(tx)
    if outgroup:
        R = t.get_common_ancestor(*outgroup)
    else:
        # Calculate the midpoint node
        R = t.get_midpoint_outgroup()

    t.set_outgroup(R)
    farthest, max_dist = t.get_farthest_leaf()

    margin = .05
    xstart = margin
    ystart = 1 - margin
    canvas = 1 - rmargin - 2 * margin
    tip = .005
    # scale the tree
    scale = canvas / max_dist

    num_leaves = len(t.get_leaf_names())
    yinterval = canvas / (num_leaves + 1)

    coords = {}
    i = 0
    for n in t.traverse("postorder"):
        dist = n.get_distance(t)
        xx = xstart + scale * dist

        if n.is_leaf():
            yy = ystart - i * yinterval
            i += 1
            ax.text(xx + tip, yy, n.name, va="center",
                    fontstyle="italic", size=8)
        else:
            children = [coords[x] for x in n.get_children()]
            children_x, children_y = zip(*children)
            min_y, max_y = min(children_y), max(children_y)
            # plot the vertical bar
            ax.plot((xx, xx), (min_y, max_y), "k-")
            # plot the horizontal bar
            for cx, cy in children:
                ax.plot((xx, cx), (cy, cy), "k-")
            yy = sum(children_y) * 1. / len(children_y)
            support = n.support
            ax.text(xx, yy, _("{0:d}".format(int(abs(support * 100)))),
                    ha="right", size=10)

        coords[n] = (xx, yy)

    # scale bar
    br = .1
    x1 = xstart + .1
    x2 = x1 + br * scale
    yy = ystart - i * yinterval
    ax.plot([x1, x1], [yy - tip, yy + tip], "k-")
    ax.plot([x2, x2], [yy - tip, yy + tip], "k-")
    ax.plot([x1, x2], [yy, yy], "k-")
    ax.text((x1 + x2) / 2, yy - tip, _("{0:g}".format(br)),
            va="top", ha="center", size=10)


def main(tx=None):
    p = OptionParser(__doc__)
    p.add_option("--outgroup",
                 help="Root the tree using the outgroup. " + \
                      "Use comma to separate multiple taxa.")
    p.add_option("--rmargin", default=.3, type="float",
                 help="Set blank rmargin to the right [default: %default]")
    opts, args, iopts = set_image_options(p, figsize="8x6")

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    outgroup = None
    if opts.outgroup:
        outgroup = opts.outgroup.split(",")
    pf = datafile.rsplit(".", 1)[0]
    if tx:
        pf = "demo"
    else:
        tx = open(datafile).read()
        logging.debug("Load tree file `{0}`.".format(datafile))

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    draw_tree(root, tx, rmargin=opts.rmargin, outgroup=outgroup)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


if __name__ == '__main__':
    t1 = """(((Os02g0681100:0.1151,Sb04g031800:0.11220)1.0:0.0537,
    (Os04g0578800:0.04318,Sb06g026210:0.04798)-1.0:0.08870)1.0:0.06985,
    ((Os03g0124100:0.08845,Sb01g048930:0.09055)1.0:0.05332,
    (Os10g0534700:0.06592,Sb01g030630:0.04824)-1.0:0.07886):0.09389);"""
    main()
