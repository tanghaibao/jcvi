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

import numpy as np
from random import sample
from itertools import groupby
from optparse import OptionParser

from jcvi.formats.bed import Bed
from jcvi.algorithms.synteny import batch_scan
from jcvi.apps.base import debug
from jcvi.graphics.base import plt, ticker, Rectangle, cm, _, \
        set_human_axis, set_image_options
debug()


def get_breaks(bed):
    # get chromosome break positions
    simple_bed = bed.simple_bed
    for seqid, ranks in groupby(simple_bed, key=lambda x: x[0]):
        ranks = list(ranks)
        # chromosome, extent of the chromosome
        yield seqid, ranks[0][1], ranks[-1][1]


def draw_box(clusters, ax, color="b"):

    for cluster in clusters:
        xrect, yrect = zip(*cluster)
        xmin, xmax, ymin, ymax = min(xrect), max(xrect), \
                                min(yrect), max(yrect)
        ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,\
                                ec=color, fc='y', alpha=.5))


def draw_cmap(ax, cmap_text, vmin, vmax, cmap=None, reverse=False):
    X = [1, 0] if reverse else [0, 1]
    Y = np.array([X, X])
    xmin, xmax = .5, .9
    ymin, ymax = .02, .04
    ax.imshow(Y, extent=(xmin, xmax, ymin, ymax), cmap=cmap)
    ax.text(xmin - .01, (ymin + ymax) * .5, _(cmap_text),
            ha="right", va="center", size=10)
    vmiddle = (vmin + vmax) * .5
    xmiddle = (xmin + xmax) * .5
    for x, v in zip((xmin, xmiddle, xmax), (vmin, vmiddle, vmax)):
        ax.text(x, ymin - .005, _("%.1f" % v), ha="center", va="top", size=10)


def dotplot(anchorfile, qbed, sbed, image_name, vmin, vmax, iopts,
        is_self=False, synteny=False, cmap_text=None):

    fp = open(anchorfile)

    qorder = qbed.order
    sorder = sbed.order

    data = []
    logging.debug("Normalize values to [%.1f, %.1f]" % (vmin, vmax))

    for row in fp:
        atoms = row.split()
        # first two columns are query and subject, and an optional third column
        if len(atoms) < 2:
            continue
        query, subject = atoms[:2]
        value = atoms[-1] if cmap_text else vmax

        try:
            value = float(value)
        except ValueError:
            value = vmin

        if value < vmin:
            value = vmin
        if value > vmax:
            value = vmax

        if query not in qorder:
            #logging.warning("ignore %s" % query)
            continue
        if subject not in sorder:
            #logging.warning("ignore %s" % subject)
            continue

        qi, q = qorder[query]
        si, s = sorder[subject]

        nv = vmax - value
        data.append((qi, si, nv))
        if is_self:  # Mirror image
            data.append((si, qi, nv))

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])  # the whole canvas
    ax = fig.add_axes([.1, .1, .8, .8])  # the dot plot

    sample_number = 5000  # only show random subset
    if len(data) > sample_number:
        data = sample(data, sample_number)

    # the data are plotted in this order, the least value are plotted
    # last for aesthetics
    data.sort(key=lambda x: -x[2])

    default_cm = cm.copper
    x, y, c = zip(*data)
    ax.scatter(x, y, c=c, s=2, lw=0, cmap=default_cm,
            vmin=vmin, vmax=vmax)

    if synteny:
        clusters = batch_scan(data, qbed, sbed)
        draw_box(clusters, ax)

    if cmap_text:
        draw_cmap(root, cmap_text, vmin, vmax, cmap=default_cm, reverse=True)

    xsize, ysize = len(qbed), len(sbed)
    logging.debug("xsize=%d ysize=%d" % (xsize, ysize))
    xlim = (0, xsize)
    ylim = (ysize, 0)  # invert the y-axis

    xchr_labels, ychr_labels = [], []
    ignore = True  # tag to mark whether to plot chr name (skip small ones)
    ignore_size_x = xsize * .005
    ignore_size_y = ysize * .005

    # plot the chromosome breaks
    for (seqid, beg, end) in get_breaks(qbed):
        ignore = abs(end - beg) < ignore_size_x
        seqid = seqid.split("_")[-1]
        try:
            seqid = int(seqid)
            seqid = "c%d" % seqid
        except:
            pass

        xchr_labels.append((seqid, (beg + end) / 2, ignore))
        ax.plot([beg, beg], ylim, "g-", lw=1)

    for (seqid, beg, end) in get_breaks(sbed):
        ignore = abs(end - beg) < ignore_size_y
        seqid = seqid.split("_")[-1]
        try:
            seqid = int(seqid)
            seqid = "c%d" % seqid
        except:
            pass

        ychr_labels.append((seqid, (beg + end) / 2, ignore))
        ax.plot(xlim, [beg, beg], "g-", lw=1)

    # plot the chromosome labels
    for label, pos, ignore in xchr_labels:
        pos = .1 + pos * .8 / xsize
        if not ignore:
            root.text(pos, .91, label,
                ha="center", va="bottom", rotation=45, color="grey")

    # remember y labels are inverted
    for label, pos, ignore in ychr_labels:
        pos = .9 - pos * .8 / ysize
        if not ignore:
            root.text(.91, pos, label,
                va="center", color="grey")

    # create a diagonal to separate mirror image for self comparison
    if is_self:
        ax.plot(xlim, (0, ysize), 'm-', alpha=.5, lw=2)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # add genome names
    to_ax_label = lambda fname: _(op.basename(fname).split(".")[0])
    gx, gy = [to_ax_label(x.filename) for x in (qbed, sbed)]
    ax.set_xlabel(gx, size=16)
    ax.set_ylabel(gy, size=16)

    # beautify the numeric axis
    for tick in ax.get_xticklines() + ax.get_yticklines():
        tick.set_visible(False)

    set_human_axis(ax)

    plt.setp(ax.get_xticklabels() + ax.get_yticklabels(),
            color='gray', size=10)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)


if __name__ == "__main__":

    p = OptionParser(__doc__)
    p.add_option("--qbed", help="Path to qbed")
    p.add_option("--sbed", help="Path to sbed")
    p.add_option("--synteny", default=False, action="store_true",
            help="Run a fast synteny scan and display blocks [default: %default]")
    p.add_option("--cmap", default="Synonymous substitutions (Ks)",
            help="Draw colormap box on the bottom-left corner "
                 "[default: `%default`]")
    p.add_option("--vmin", dest="vmin", type="float", default=0,
            help="Minimum value in the colormap [default: %default]")
    p.add_option("--vmax", dest="vmax", type="float", default=1,
            help="Maximum value in the colormap [default: %default]")
    opts, args, iopts = set_image_options(p, figsize="8x8", dpi=90)

    qbed, sbed = opts.qbed, opts.sbed
    if not (len(args) == 1 and qbed and sbed):
        sys.exit(p.print_help())

    is_self = (qbed == sbed)

    qbed = Bed(qbed)
    sbed = Bed(sbed)
    synteny = opts.synteny
    vmin, vmax = opts.vmin, opts.vmax
    cmap_text = opts.cmap

    anchorfile = args[0]

    image_name = op.splitext(anchorfile)[0] + "." + opts.format
    dotplot(anchorfile, qbed, sbed, image_name, vmin, vmax, iopts,
            is_self=is_self, synteny=synteny, cmap_text=cmap_text)
