#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [anchorfile|ksfile] --qbed query.bed --sbed subject.bed

visualize the anchorfile in a dotplot. anchorfile contains two columns
indicating gene pairs, followed by an optional column (e.g. Ks value).

The option --colormap specifies the block color to highlight certain blocks in
a file.  Block ids are 1-based (non-digit chars will be removed). For example, below
requests that block 1 is class 'sigma' and block 2 is class 'tau'.

1    sigma
2    tau
3    tau

These classes will be mapped to auto-assigned colors and figure legend added to
the bottom of the figure.

*Important*

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

from jcvi.compara.synteny import AnchorFile, batch_scan, check_beds, get_orientation
from jcvi.utils.cbook import seqid_parse, thousands
from jcvi.apps.base import OptionParser, need_update
from jcvi.graphics.base import (
    plt,
    Rectangle,
    set_human_axis,
    savefig,
    draw_cmap,
    TextHandler,
    latex,
    markup,
    normalize_axes,
    set1,
)


class Palette(dict):
    def __init__(self, palettedict=None, palettefile=None):
        """ Instantiate a palette to map from block_id to color

        Args:
            palettedict (Dict, optional): Get the mapping from a dict. Defaults to None.
            palettefile (str, optional): Get the mapping from a two-column file. Defaults to None.
        """
        if palettedict is not None:
            self.update(palettedict)
        if palettefile is None:
            return

        pal = "rbcygmk"

        fp = open(palettefile)
        for row in fp:
            a, b = row.split()
            a = "".join(x for x in a if x in string.digits)
            a = int(a)
            self[a] = b

        self.categories = sorted(set(self.values()))
        self.colors = dict(zip(self.categories, pal))

        logging.debug(
            "Color info ({0} categories) imported for {1} blocks.".format(
                len(self.colors), len(self)
            )
        )
        logging.debug(str(self.colors))

        for k, v in self.items():  # Update from categories to colors
            self[k] = self.colors[v]

    @classmethod
    def from_block_orientation(
        cls, anchorfile, qbed, sbed, forward_color="#e7298a", reverse_color="#3690c0"
    ):
        """ Generate a palette which contains mapping from block_id (1-based) to colors.

        Args:
            anchorfile (str): Path to the .anchors file
            qbed (BedFile): Query BED
            sbed (BedFile): Subject BED
            forward_color (str, optional): Color of forward block. Defaults to "#e7298a".
            reverse_color (str, optional): Color of reverse block. Defaults to "#3690c0".
        """
        ac = AnchorFile(anchorfile)
        blocks = ac.blocks
        palette = {}
        qorder = qbed.order
        sorder = sbed.order

        for i, block in enumerate(blocks):
            block_id = i + 1

            a, b, _ = zip(*block)
            a = [qorder[x] for x in a]
            b = [sorder[x] for x in b]
            ia, _ = zip(*a)
            ib, _ = zip(*b)

            orientation = get_orientation(ia, ib)
            palette[block_id] = reverse_color if orientation == "-" else forward_color
        return cls(palettedict=palette)


def draw_box(clusters, ax, color="b"):

    for cluster in clusters:
        xrect, yrect = zip(*cluster)
        xmin, xmax, ymin, ymax = min(xrect), max(xrect), min(yrect), max(yrect)
        ax.add_patch(
            Rectangle(
                (xmin, ymin), xmax - xmin, ymax - ymin, ec=color, fc="y", alpha=0.5
            )
        )


def plot_breaks_and_labels(
    fig,
    root,
    ax,
    gx,
    gy,
    xsize,
    ysize,
    qbreaks,
    sbreaks,
    sep=True,
    chrlw=0.1,
    sepcolor="g",
    minfont=5,
    stdpf=True,
    chpf=True,
):
    xlim = (0, xsize)
    ylim = (ysize, 0)  # invert the y-axis

    # Tag to mark whether to plot chr name (skip small ones)
    xchr_labels, ychr_labels = [], []
    th = TextHandler(fig)

    # plot the chromosome breaks
    for (seqid, beg, end) in qbreaks:
        xsize_ratio = abs(end - beg) * 0.8 / xsize
        fontsize = th.select_fontsize(xsize_ratio)
        if chpf:
            seqid = "".join(seqid_parse(seqid, stdpf=stdpf)[:2])

        xchr_labels.append((seqid, (beg + end) / 2, fontsize))
        if sep:
            ax.plot([beg, beg], ylim, "-", lw=chrlw, color=sepcolor)

    for (seqid, beg, end) in sbreaks:
        ysize_ratio = abs(end - beg) * 0.8 / ysize
        fontsize = th.select_fontsize(ysize_ratio)
        if chpf:
            seqid = "".join(seqid_parse(seqid, stdpf=stdpf)[:2])

        ychr_labels.append((seqid, (beg + end) / 2, fontsize))
        if sep:
            ax.plot(xlim, [beg, beg], "-", lw=chrlw, color=sepcolor)

    # plot the chromosome labels
    for label, pos, fontsize in xchr_labels:
        pos = 0.1 + pos * 0.8 / xsize
        if fontsize >= minfont:
            root.text(
                pos,
                0.91,
                latex(label),
                size=fontsize,
                ha="center",
                va="bottom",
                rotation=45,
                color="grey",
            )

    # remember y labels are inverted
    for label, pos, fontsize in ychr_labels:
        pos = 0.9 - pos * 0.8 / ysize
        if fontsize >= minfont:
            root.text(0.91, pos, latex(label), size=fontsize, va="center", color="grey")

    # Plot the frame
    ax.plot(xlim, [0, 0], "-", lw=chrlw, color=sepcolor)
    ax.plot(xlim, [ysize, ysize], "-", lw=chrlw, color=sepcolor)
    ax.plot([0, 0], ylim, "-", lw=chrlw, color=sepcolor)
    ax.plot([xsize, xsize], ylim, "-", lw=chrlw, color=sepcolor)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.set_xlabel(gx, size=16)
    ax.set_ylabel(gy, size=16)

    # beautify the numeric axis
    for tick in ax.get_xticklines() + ax.get_yticklines():
        tick.set_visible(False)

    set_human_axis(ax)

    plt.setp(ax.get_xticklabels() + ax.get_yticklabels(), color="gray", size=10)

    return xlim, ylim


def downsample(data, sample_number=10000):
    npairs = len(data)
    # Only show random subset
    if npairs > sample_number:
        logging.debug(
            "Showing a random subset of {0} data points (total {1}) "
            "for clarity.".format(sample_number, npairs)
        )
        data = sample(data, sample_number)
    return npairs


def dotplot(
    anchorfile,
    qbed,
    sbed,
    fig,
    root,
    ax,
    vmin=0,
    vmax=1,
    is_self=False,
    synteny=False,
    cmap_text=None,
    cmap="copper",
    genomenames=None,
    sample_number=10000,
    minfont=5,
    palette=None,
    chrlw=0.1,
    title=None,
    sep=True,
    sepcolor="g",
    stdpf=True,
    chpf=True,
):

    fp = open(anchorfile)
    # add genome names
    if genomenames:
        gx, gy = genomenames.split("_")
    else:
        to_ax_label = lambda fname: op.basename(fname).split(".")[0]
        gx, gy = [to_ax_label(x.filename) for x in (qbed, sbed)]

    # Stylize the axis labels
    gx, gy = markup(gx), markup(gy)

    qorder = qbed.order
    sorder = sbed.order

    data = []
    if cmap_text:
        logging.debug("Capping values within [{0:.1f}, {1:.1f}]".format(vmin, vmax))

    block_id = 0
    for row in fp:
        atoms = row.split()
        if row[0] == "#":
            block_id += 1
            block_color = palette.get(block_id, "k") if palette else None
            continue

        # first two columns are query and subject, and an optional third column
        if len(atoms) < 2:
            continue

        query, subject = atoms[:2]
        value = atoms[-1]

        if cmap_text:
            try:
                value = float(value)
            except ValueError:
                value = vmax

            if value < vmin:
                continue
            if value > vmax:
                continue
        else:
            value = 0

        if query not in qorder:
            continue
        if subject not in sorder:
            continue

        qi, q = qorder[query]
        si, s = sorder[subject]

        nv = block_color or value
        data.append((qi, si, nv))
        if is_self:  # Mirror image
            data.append((si, qi, nv))

    npairs = downsample(data, sample_number=sample_number)
    x, y, c = zip(*data)

    if palette:
        ax.scatter(x, y, c=c, edgecolors="none", s=2, lw=0)
    else:
        ax.scatter(
            x, y, c=c, edgecolors="none", s=2, lw=0, cmap=cmap, vmin=vmin, vmax=vmax
        )

    if synteny:
        clusters = batch_scan(data, qbed, sbed)
        draw_box(clusters, ax)

    if cmap_text:
        draw_cmap(root, cmap_text, vmin, vmax, cmap=cmap)

    xsize, ysize = len(qbed), len(sbed)
    logging.debug("xsize=%d ysize=%d" % (xsize, ysize))
    qbreaks = qbed.get_breaks()
    sbreaks = sbed.get_breaks()
    xlim, ylim = plot_breaks_and_labels(
        fig,
        root,
        ax,
        gx,
        gy,
        xsize,
        ysize,
        qbreaks,
        sbreaks,
        sep=sep,
        chrlw=chrlw,
        sepcolor=sepcolor,
        minfont=minfont,
        stdpf=stdpf,
        chpf=chpf,
    )

    # create a diagonal to separate mirror image for self comparison
    if is_self:
        ax.plot(xlim, (0, ysize), "m-", alpha=0.5, lw=2)

    if palette and hasattr(
        palette, "colors"
    ):  # bottom-left has the palette, if available
        colors = palette.colors
        xstart, ystart = 0.1, 0.05
        for category, c in sorted(colors.items()):
            root.add_patch(Rectangle((xstart, ystart), 0.03, 0.02, lw=0, fc=c))
            root.text(xstart + 0.04, ystart, category, color=c)
            xstart += 0.1

    if title is None:
        title = "Inter-genomic comparison: {0} vs {1}".format(gx, gy)
        if is_self:
            title = "Intra-genomic comparison within {0}".format(gx)
            npairs /= 2
        title += " ({0} gene pairs)".format(thousands(npairs))
    root.set_title(title, x=0.5, y=0.96, color="k")
    if title:
        logging.debug("Dot plot title: {}".format(title))
    normalize_axes(root)


def subset_bed(bed, seqids):
    from copy import deepcopy

    newbed = deepcopy(bed)
    del newbed[:]
    for b in bed:
        if b.seqid not in seqids:
            continue
        newbed.append(b)
    return newbed


def dotplot_main(args):
    p = OptionParser(__doc__)
    p.set_beds()
    p.add_option(
        "--synteny",
        default=False,
        action="store_true",
        help="Run a fast synteny scan and display blocks",
    )
    p.add_option("--cmaptext", help="Draw colormap box on the bottom-left corner")
    p.add_option(
        "--vmin",
        dest="vmin",
        type="float",
        default=0,
        help="Minimum value in the colormap",
    )
    p.add_option(
        "--vmax",
        dest="vmax",
        type="float",
        default=2,
        help="Maximum value in the colormap",
    )
    p.add_option(
        "--nmax",
        dest="sample_number",
        type="int",
        default=10000,
        help="Maximum number of data points to plot",
    )
    p.add_option(
        "--minfont",
        type="int",
        default=4,
        help="Do not render labels with size smaller than",
    )
    p.add_option("--colormap", help="Two column file, block id to color mapping")
    p.add_option(
        "--colororientation",
        action="store_true",
        default=False,
        help="Color the blocks based on orientation, similar to mummerplot",
    )
    p.add_option(
        "--nosort",
        default=False,
        action="store_true",
        help="Do not sort the seqids along the axes",
    )
    p.add_option(
        "--nosep", default=False, action="store_true", help="Do not add contig lines"
    )
    p.add_option("--title", help="Title of the dot plot")
    p.set_dotplot_opts()
    p.set_outfile(outfile=None)
    opts, args, iopts = p.set_image_options(
        args, figsize="9x9", style="dark", dpi=90, cmap="copper"
    )

    if len(args) != 1:
        sys.exit(not p.print_help())

    (anchorfile,) = args
    qbed, sbed, qorder, sorder, is_self = check_beds(
        anchorfile, p, opts, sorted=(not opts.nosort)
    )

    palette = opts.colormap
    if palette:
        palette = Palette(palettefile=palette)
    elif opts.colororientation:
        palette = Palette.from_block_orientation(anchorfile, qbed, sbed)

    cmaptext = opts.cmaptext
    if anchorfile.endswith(".ks"):
        from jcvi.apps.ks import KsFile

        logging.debug("Anchors contain Ks values")
        cmaptext = cmaptext or "*Ks* values"
        anchorksfile = anchorfile + ".anchors"
        if need_update(anchorfile, anchorksfile):
            ksfile = KsFile(anchorfile)
            ksfile.print_to_anchors(anchorksfile)
        anchorfile = anchorksfile

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

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])  # the whole canvas
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # the dot plot

    dotplot(
        anchorfile,
        qbed,
        sbed,
        fig,
        root,
        ax,
        vmin=opts.vmin,
        vmax=opts.vmax,
        is_self=is_self,
        synteny=opts.synteny,
        cmap_text=opts.cmaptext,
        cmap=iopts.cmap,
        genomenames=opts.genomenames,
        sample_number=opts.sample_number,
        minfont=opts.minfont,
        palette=palette,
        sep=(not opts.nosep),
        sepcolor=set1[int(opts.theme)],
        title=opts.title,
        stdpf=(not opts.nostdpf),
        chpf=(not opts.nochpf),
    )

    image_name = opts.outfile or (op.splitext(anchorfile)[0] + "." + opts.format)
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)
    fig.clear()


if __name__ == "__main__":
    dotplot_main(sys.argv[1:])
