#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog anchorfile --qbed query.bed --sbed subject.bed

visualize the anchorfile in a dotplot. anchorfile contains two columns
indicating gene pairs, followed by an optional column (e.g. Ks value)
"""

import os.path as op
import sys
import logging
from itertools import groupby
from optparse import OptionParser

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

from jcvi.formats.bed import Bed
from jcvi.algorithms.synteny import batch_scan
from jcvi.apps.base import debug
debug()


def get_breaks(bed):
    # get chromosome break positions
    simple_bed = bed.simple_bed
    for seqid, ranks in groupby(simple_bed, key=lambda x:x[0]):
        ranks = list(ranks)
        # chromosome, extent of the chromosome
        yield seqid, ranks[0][1], ranks[-1][1]


def draw_box(clusters, ax, color="b"):

    for cluster in clusters:
        xrect, yrect = zip(*cluster)
        xmin, xmax, ymin, ymax = min(xrect), max(xrect), \
                                min(yrect), max(yrect)
        ax.add_patch(Rectangle((xmin, ymin), xmax-xmin, ymax-ymin,\
                                ec=color, fc='y', alpha=.5))
                                

def dotplot(anchorfile, qbed, sbed, image_name, is_self=False, synteny=False):

    fp = open(anchorfile)

    qorder = qbed.order
    sorder = sbed.order

    data = []
    for row in fp:
        atoms = row.split()
        if len(atoms) < 2: continue
        query, subject = row.split()[:2]
        #value = float(value)
        if query not in qorder: 
            #logging.warning("ignore %s" % query)
            continue
        if subject not in sorder: 
            #logging.warning("ignore %s" % subject)
            continue

        qi, q = qorder[query]
        si, s = sorder[subject]
        data.append((qi, si))

    fig = plt.figure(1,(8,8))
    root = fig.add_axes([0,0,1,1]) # the whole canvas
    ax = fig.add_axes([.1,.1,.8,.8]) # the dot plot

    x, y = zip(*data)
    ax.scatter(x, y, c='k', s=.05, lw=0, alpha=.9)

    if synteny:
        clusters = batch_scan(data, qbed, sbed)
        draw_box(clusters, ax)

    xsize, ysize = len(qbed), len(sbed)
    xlim = (0, xsize)
    ylim = (ysize, 0) # invert the y-axis

    xchr_labels, ychr_labels = [], []
    ignore = True # tag to mark whether to plot chromosome name (skip small ones)
    ignore_size = 100
    # plot the chromosome breaks
    for (seqid, beg, end) in get_breaks(qbed):
        ignore = abs(end-beg) < ignore_size 
        xchr_labels.append((seqid, (beg + end)/2, ignore))
        ax.plot([beg, beg], ylim, "g-", alpha=.5)

    for (seqid, beg, end) in get_breaks(sbed):
        ignore = abs(end-beg) < ignore_size 
        ychr_labels.append((seqid, (beg + end)/2, ignore))
        ax.plot(xlim, [beg, beg], "g-", alpha=.5)

    # plot the chromosome labels
    for label, pos, ignore in xchr_labels:
        pos = .1 + pos * .8/ xsize
        if not ignore:
            root.text(pos, .91, r"%s" % label, color="b",
                size=9, alpha=.5, rotation=45)

    # remember y labels are inverted
    for label, pos, ignore in ychr_labels:
        pos = .9 - pos * .8/ ysize
        if not ignore:
            root.text(.91, pos, r"%s" % label, color="b",
                size=9, alpha=.5, ha="left", va="center")

    # create a diagonal to separate mirror image for self comparison
    if is_self:
        ax.plot(xlim, ylim, 'm-', alpha=.5, lw=2)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # i always like the latex font
    _ = lambda x: r"$\mathsf{%s}$" % x.replace("_", " ").replace(" ", r"\ ")
    to_ax_label = lambda fname: _(op.basename(fname).split(".")[0])

    # add genome names
    ax.set_xlabel(to_ax_label(qbed.filename))
    ax.set_ylabel(to_ax_label(sbed.filename))

    # beautify the numeric axis
    for tick in ax.get_xticklines() + ax.get_yticklines():
        tick.set_visible(False) 

    formatter = ticker.FuncFormatter(lambda x, pos: r"$%d$" % x)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

    plt.setp(ax.get_xticklabels() + ax.get_yticklabels(), color='gray', size=10)

    root.set_axis_off()
    logging.debug("print image to %s" % image_name)
    plt.savefig(image_name, dpi=1000)


if __name__ == "__main__":

    parser = OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed",
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed",
            help="path to sbed")
    parser.add_option("--synteny", dest="synteny", 
            default=False, action="store_true",
            help="run a fast synteny scan and display synteny blocks")
    parser.add_option("--format", dest="format", default="png",
            help="generate image of format (png, pdf, ps, eps, svg, etc.)"
            "[default: %default]")

    opts, args = parser.parse_args()

    qbed, sbed = opts.qbed, opts.sbed
    if not (len(args) == 1 and qbed and sbed):
        sys.exit(parser.print_help())

    is_self = False
    if qbed==sbed:
        print >>sys.stderr, "Looks like this is self-self comparison"
        is_self = True

    qbed = Bed(qbed)
    sbed = Bed(sbed)
    synteny = opts.synteny

    anchorfile = args[0]

    image_name = op.splitext(anchorfile)[0] + "." + opts.format
    dotplot(anchorfile, qbed, sbed, image_name, is_self=is_self, synteny=synteny)

