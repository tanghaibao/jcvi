#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog anchorfile --qbed query.bed --sbed subject.bed

visualize the anchorfile in a dotplot. anchorfile contains two columns
indicating gene pairs, followed by an optional column (e.g. Ks value).

The option --colormap specifies the block color to highlight certain blocks in
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

from random import sample

from jcvi.compara.synteny import AnchorFile, batch_scan, check_beds
from jcvi.utils.cbook import seqid_parse, thousands
from jcvi.apps.base import OptionParser
from jcvi.graphics.base import plt, Rectangle, cm, set_human_axis, savefig, \
            draw_cmap, TextHandler, latex


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


def dotplot_main(anchorfile, qbed, sbed, image_name, iopts, vmin=0, vmax=1,
        is_self=False, synteny=False, cmap_text=None, genomenames=None,
        sample_number=10000, minfont=5, palette=None, chrlw=.01, title=None):

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])  # the whole canvas
    ax = fig.add_axes([.1, .1, .8, .8])  # the dot plot

    dotplot(anchorfile, qbed, sbed, fig, root, ax, vmin=vmin, vmax=vmax,
        is_self=is_self, synteny=synteny, cmap_text=cmap_text,
        genomenames=genomenames, sample_number=sample_number,
        minfont=minfont, palette=palette, chrlw=chrlw, title=title)

    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def dotplot(anchorfile, qbed, sbed, fig, root, ax, vmin=0, vmax=1,
        is_self=False, synteny=False, cmap_text=None, genomenames=None,
        sample_number=10000, minfont=5, palette=None, chrlw=.01, title=None,
        sepcolor="gainsboro"):

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
            continue
        if subject not in sorder:
            continue

        qi, q = qorder[query]
        si, s = sorder[subject]

        nv = vmax - value if block_color is None else block_color
        data.append((qi, si, nv))
        if is_self:  # Mirror image
            data.append((si, qi, nv))

    npairs = len(data)
    # Only show random subset
    if npairs > sample_number:
        logging.debug("Showing a random subset of {0} data points (total {1}) " \
                      "for clarity.".format(sample_number, npairs))
        data = sample(data, sample_number)

    # the data are plotted in this order, the least value are plotted
    # last for aesthetics
    if not palette:
        data.sort(key=lambda x: -x[2])

    default_cm = cm.copper
    x, y, c = zip(*data)

    if palette:
        ax.scatter(x, y, c=c, edgecolors="none", s=2, lw=0)

    else:
        ax.scatter(x, y, c=c, edgecolors="none", s=2, lw=0, cmap=default_cm,
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

    # Tag to mark whether to plot chr name (skip small ones)
    xchr_labels, ychr_labels = [], []
    th = TextHandler(fig)

    # plot the chromosome breaks
    for (seqid, beg, end) in qbed.get_breaks():
        xsize_ratio = abs(end - beg) * .8 / xsize
        fontsize = th.select_fontsize(xsize_ratio)
        seqid = "".join(seqid_parse(seqid)[:2])

        xchr_labels.append((seqid, (beg + end) / 2, fontsize))
        ax.plot([beg, beg], ylim, "-", lw=chrlw, color=sepcolor)

    for (seqid, beg, end) in sbed.get_breaks():
        ysize_ratio = abs(end - beg) * .8 / ysize
        fontsize = th.select_fontsize(ysize_ratio)
        seqid = "".join(seqid_parse(seqid)[:2])

        ychr_labels.append((seqid, (beg + end) / 2, fontsize))
        ax.plot(xlim, [beg, beg], "-", lw=chrlw, color=sepcolor)

    # plot the chromosome labels
    for label, pos, fontsize in xchr_labels:
        pos = .1 + pos * .8 / xsize
        if fontsize >= minfont:
            root.text(pos, .91, latex(label), size=fontsize,
                ha="center", va="bottom", rotation=45, color="grey")

    # remember y labels are inverted
    for label, pos, fontsize in ychr_labels:
        pos = .9 - pos * .8 / ysize
        if fontsize >= minfont:
            root.text(.91, pos, latex(label), size=fontsize,
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
        to_ax_label = lambda fname: op.basename(fname).split(".")[0]
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

    if not title:
        title = "Inter-genomic comparison: {0} vs {1}".format(gx, gy)
        if is_self:
            title = "Intra-genomic comparison within {0}".format(gx)
            npairs /= 2
        title += " ({0} gene pairs)".format(thousands(npairs))
    root.set_title(title, x=.5, y=.96, color="k")
    logging.debug(title)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()


def subset_bed(bed, seqids):
    from copy import deepcopy

    newbed = deepcopy(bed)
    del newbed[:]
    for b in bed:
        if b.seqid not in seqids:
            continue
        newbed.append(b)
    return newbed


if __name__ == "__main__":

    p = OptionParser(__doc__)
    p.set_beds()
    p.add_option("--synteny", default=False, action="store_true",
            help="Run a fast synteny scan and display blocks [default: %default]")
    p.add_option("--cmaptext", help="Draw colormap box on the bottom-left corner")
    p.add_option("--vmin", dest="vmin", type="float", default=0,
            help="Minimum value in the colormap [default: %default]")
    p.add_option("--vmax", dest="vmax", type="float", default=1,
            help="Maximum value in the colormap [default: %default]")
    p.add_option("--genomenames", type="string", default=None,
            help="genome names for labeling axes in the form of qname_sname, " \
            "eg. \"Vitis vinifera_Oryza sativa\"")
    p.add_option("--nmax", dest="sample_number", type="int", default=10000,
            help="Maximum number of data points to plot [default: %default]")
    p.add_option("--minfont", type="int", default=4,
            help="Do not render labels with size smaller than")
    p.add_option("--colormap",
            help="Two column file, block id to color mapping [default: %default]")
    p.add_option("--skipempty", default=False, action="store_true",
            help="Skip seqids that do not have matches")
    opts, args, iopts = p.set_image_options(sys.argv[1:], figsize="8x8",
                                            style="dark", dpi=90)

    if len(args) != 1:
        sys.exit(not p.print_help())

    palette = opts.colormap
    if palette:
        palette = Palette(palette)

    anchorfile, = args
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)

    if opts.skipempty:
        ac = AnchorFile(anchorfile)
        if is_self:
            qseqids = sseqids = set()
        else:
            qseqids, sseqids = set(), set()

        for pair in ac.iter_pairs():
            q, s = pair[:2]
            qi, q = qorder[q]
            si, s = sorder[s]
            qseqids.add(q.seqid)
            sseqids.add(s.seqid)

        if is_self:
            qbed = sbed = subset_bed(qbed, qseqids)
        else:
            qbed = subset_bed(qbed, qseqids)
            sbed = subset_bed(sbed, sseqids)

    image_name = op.splitext(anchorfile)[0] + "." + opts.format
    dotplot_main(anchorfile, qbed, sbed, image_name, iopts,
            vmin=opts.vmin, vmax=opts.vmax, is_self=is_self,
            synteny=opts.synteny, cmap_text=opts.cmaptext,
            genomenames=opts.genomenames, sample_number=opts.sample_number,
            minfont=opts.minfont, palette=palette)
