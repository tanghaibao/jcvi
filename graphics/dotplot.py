#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog anchorfile --qbed query.bed --sbed subject.bed

visualize the anchorfile in a dotplot. anchorfile contains two columns
indicating gene pairs, followed by an optional column (e.g. Ks value).

The option --palette specifies the block color to highlight certain blocks in
a file.  Block ids are 1-based (non-digit chars will be removed). For example, below
requests that the 7th blocks to be colored red.

rice-sigma07    sigma
rice-sigma10    tau

Before running this script it is recommended to check/install
TeX Live (http://www.tug.org/texlive/) and
Ghostscript (http://www.ghostscript.com/)
see more here: http://matplotlib.sourceforge.net/users/usetex.html
"""

import os.path as op
import sys
import logging
import string

import numpy as np
from random import sample
from itertools import groupby
from optparse import OptionParser

from jcvi.formats.bed import Bed
from jcvi.algorithms.synteny import batch_scan, add_beds, check_beds
from jcvi.apps.base import debug
from jcvi.graphics.base import plt, ticker, Rectangle, cm, _, \
        set_human_axis, set_image_options, savefig
debug()


class Palette (dict):

    def __init__(self, palettefile):

        pal = "rbcygmk"

        fp = open(palettefile)
        for row in fp:
            a, b = row.split()
            a = "".join(x for x in a if x in string.digits)
            a = int(a)
            self[a] = b

        self.categories = sorted(set(self.values()))
        self.colors = dict(zip(self.categories, pal))

        logging.debug("Color info ({0} categories) imported for {1} blocks.".\
                        format(len(self.colors), len(self)))
        logging.debug(str(self.colors))

        for k, v in self.items():  # Update from categories to colors
            self[k] = self.colors[v]


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
        is_self=False, synteny=False, cmap_text=None, genomenames=None,
        sample_number=10000, ignore=.005, palette=None, chrlw=1):

    fp = open(anchorfile)

    qorder = qbed.order
    sorder = sbed.order

    data = []
    if cmap_text:
        logging.debug("Normalize values to [%.1f, %.1f]" % (vmin, vmax))

    block_id = 0
    for row in fp:
        atoms = row.split()
        block_color = None
        if row[0] == "#":
            block_id += 1
            if palette:
                block_color = palette.get(block_id, "k")
            continue

        # first two columns are query and subject, and an optional third column
        if len(atoms) < 2:
            continue

        query, subject = atoms[:2]
        value = atoms[-1]

        try:
            value = float(value)
        except ValueError:
            value = vmax

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

        nv = vmax - value if block_color is None else block_color
        data.append((qi, si, nv))
        if is_self:  # Mirror image
            data.append((si, qi, nv))

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])  # the whole canvas
    ax = fig.add_axes([.1, .1, .8, .8])  # the dot plot

    # only show random subset, default to sample_number = 5000
    if len(data) > sample_number:
        logging.debug("Showing a random subset of %s data points (total %s) " \
                      "for clarity." % (sample_number, len(data)))
        data = sample(data, sample_number)

    # the data are plotted in this order, the least value are plotted
    # last for aesthetics
    if not palette:
        data.sort(key=lambda x: -x[2])

    default_cm = cm.copper
    x, y, c = zip(*data)

    if palette:
        ax.scatter(x, y, c=c, s=2, lw=0)

    else:
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
    # Tag to mark whether to plot chr name (skip small ones)
    ignore_size_x = ignore_size_y = 0
    if ignore:
        ignore_size_x = xsize * ignore
        ignore_size_y = ysize * ignore

    # plot the chromosome breaks
    for (seqid, beg, end) in qbed.get_breaks():
        ignore = abs(end - beg) < ignore_size_x
        seqid = seqid.split("_")[-1]
        try:
            seqid = int(seqid)
            seqid = "c%d" % seqid
        except:
            pass

        xchr_labels.append((seqid, (beg + end) / 2, ignore))
        ax.plot([beg, beg], ylim, "g-", lw=chrlw)

    for (seqid, beg, end) in sbed.get_breaks():
        ignore = abs(end - beg) < ignore_size_y
        seqid = seqid.split("_")[-1]
        try:
            seqid = int(seqid)
            seqid = "c%d" % seqid
        except:
            pass

        ychr_labels.append((seqid, (beg + end) / 2, ignore))
        ax.plot(xlim, [beg, beg], "g-", lw=chrlw)

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
    if genomenames:
        gx, gy = genomenames.split("_")
    else:
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

    if palette:  # bottom-left has the palette, if available
        colors = palette.colors
        xstart, ystart = .1, .05
        for category, c in sorted(colors.items()):
            root.add_patch(Rectangle((xstart, ystart), .03, .02, lw=0, fc=c))
            root.text(xstart + .04, ystart, category, color=c)
            xstart += .1

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":

    p = OptionParser(__doc__)
    add_beds(p)
    p.add_option("--synteny", default=False, action="store_true",
            help="Run a fast synteny scan and display blocks [default: %default]")
    p.add_option("--cmap", help="Draw colormap box on the bottom-left corner "
                 "[default: `%default`]")
    p.add_option("--vmin", dest="vmin", type="float", default=0,
            help="Minimum value in the colormap [default: %default]")
    p.add_option("--vmax", dest="vmax", type="float", default=1,
            help="Maximum value in the colormap [default: %default]")
    p.add_option("--genomenames", type="string", default=None,
            help="genome names for labeling axes in the form of qname_sname, " \
            "eg. \"Vitis vinifera_Oryza sativa\"")
    p.add_option("--nmax", dest="sample_number", type="int", default=10000,
            help="Maximum number of data points to plot [default: %default]")
    p.add_option("--ignore", type="float", default=.005,
            help="Do not render labels for chr less than portion of genome [default: %default]")
    p.add_option("--palette",
            help="Two column file, block id to color mapping [default: %default]")
    opts, args, iopts = set_image_options(p, sys.argv[1:], figsize="8x8", dpi=90)

    if len(args) != 1:
        sys.exit(not p.print_help())

    synteny = opts.synteny
    vmin, vmax = opts.vmin, opts.vmax
    cmap_text = opts.cmap
    genomenames = opts.genomenames
    sample_number = opts.sample_number
    palette = opts.palette
    if palette:
        palette = Palette(palette)

    anchorfile, = args
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)

    image_name = op.splitext(anchorfile)[0] + "." + opts.format
    dotplot(anchorfile, qbed, sbed, image_name, vmin, vmax, iopts, \
            is_self=is_self, synteny=synteny, cmap_text=cmap_text, \
            genomenames=genomenames, sample_number=sample_number,
            ignore=opts.ignore, palette=palette)
