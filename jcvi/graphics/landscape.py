#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Create chromosome landscape plots that are similar to the ones used in soybean
and sorghum paper.
"""


import os.path as op
import sys
import logging

import numpy as np

from collections import defaultdict

from jcvi.formats.sizes import Sizes
from jcvi.formats.base import LineFile, DictFile
from jcvi.formats.bed import Bed, bins
from jcvi.algorithms.matrix import moving_sum
from jcvi.graphics.base import (
    plt,
    Rectangle,
    CirclePolygon,
    savefig,
    ticker,
    human_readable_base,
    latex,
    set_human_axis,
)
from jcvi.utils.cbook import human_size, autoscale
from jcvi.apps.base import OptionParser, ActionDispatcher


# Colors picked from Schmutz soybean genome paper using ColorPic
palette = ["#ACABD5", "#DBF0F5", "#3EA77A", "#FBF5AB", "#C162A6"] + list("rgbymck")
gray = "#CCCCCB"
Registration = {
    "Gypsy": "LTR-RT/Gypsy",
    "Copia": "LTR-RT/Copia",
    "hAT": "DNA-TE/hAT",
    "Helitron": "DNA-TE/Helitron",
    "Tourist": "DNA-TE/Tourist",
    "Introns": "Genes (introns)",
    "Exons": "Genes (exons)",
}


class BinLine(object):
    def __init__(self, row):
        args = row.split()
        self.chr = args[0]
        self.len = float(args[1])
        self.binlen = int(args[2])

    def __str__(self):
        return "\t".join(str(x) for x in (self.chr, self.len, self.binlen))

    def subtract(self, o):
        self.binlen -= o.len


class BinFile(LineFile):
    def __init__(self, filename):
        super(BinFile, self).__init__(filename)
        self.mapping = defaultdict(list)

        fp = open(filename)
        for row in fp:
            b = BinLine(row)
            self.append(b)
            chr, len, binlen = b.chr, b.len, b.binlen
            self.mapping[chr].append((len, binlen))
        fp.close()


def main():

    actions = (
        ("stack", "create landscape plot with genic/te composition"),
        ("heatmap", "similar to stack but adding heatmap"),
        ("composite", "combine line plots, feature bars and alt-bars"),
        ("multilineplot", "combine multiple line plots in one vertical stack"),
        ("depth", "show per chromosome depth plot across genome"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def draw_depth(root, ax, bed, colormap=None, defaultcolor="k"):
    """ Draw depth plot on the given axes, using data from bed

    Args:
        root (matplotlib.Axes): Canvas axes
        ax (matplotlib.Axes): Axes to plot data on
        bed (Bed): Bed data from mosdepth
        colormap (Dict): seqid => color
    """
    if colormap is None:
        colormap = {}
    sizes = bed.max_bp_in_chr
    seqids = colormap.keys() if colormap else sizes.keys()
    starts = {}
    label_positions = []
    start = 0
    for seqid in seqids:
        starts[seqid] = start
        end = start + sizes[seqid]
        start = end
        label_positions.append((seqid, (start + end) / 2))
    xsize = end

    # Extract plotting data
    data = []
    for b in bed:
        seqid = b.seqid
        if seqid not in starts:
            continue
        # chr01A  2000000 3000000 113.00
        x = starts[seqid] + (b.start + b.end) / 2
        y = float(b.accn)
        c = colormap.get(seqid, defaultcolor)
        data.append((x, y, c))

    x, y, c = zip(*data)
    ax.scatter(
        x, y, c=c, edgecolors="none", s=2, lw=0,
    )
    logging.debug("Obtained {} data points with depth data".format(len(data)))

    # beautify the numeric axis
    for tick in ax.get_xticklines() + ax.get_yticklines():
        tick.set_visible(False)

    for seqid, position in label_positions:
        xpos = 0.1 + position * 0.8 / xsize
        c = colormap.get(seqid, defaultcolor)
        root.text(xpos, 0.05, seqid, color=c, ha="right", va="center", rotation=20)

    ax.set_xlim(0, xsize)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Depth")

    set_human_axis(ax)
    plt.setp(ax.get_xticklabels() + ax.get_yticklabels(), color="gray", size=10)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()


def depth(args):
    """
    %prog depth sample.regions.bed.gz

    Plot the mosdepth regions BED file. We recommend to generate this BED file
    by (please adjust the --by parameter to your required resolution):

    $ mosdepth --no-per-base --use-median --fast-mode --by 1000000 sample.wgs
    sample.bam

    Use --colormap to specify a colormap between seqid and desired color. For
    example:

    chr1    r
    chr2    g
    ...

    Only seqids that are in the colormap will be plotted, in the order that's
    given in the file. When --colormap is not set, every seqid will be drawn in black.
    """
    p = OptionParser(depth.__doc__)
    p.add_option(
        "--colormap", help="Two-column tab-separated mappings between seqid and color"
    )
    opts, args, iopts = p.set_image_options(args, style="dark", figsize="12x4")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    bed = Bed(bedfile)

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    ax = fig.add_axes([0.1, 0.2, 0.8, 0.7])

    colormap = DictFile(opts.colormap) if opts.colormap else {}
    draw_depth(root, ax, bed, colormap=colormap)

    pf = bedfile.split(".", 1)[0]
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def add_window_options(p):
    p.add_option("--window", default=500000, type="int", help="Size of window")
    p.add_option("--shift", default=100000, type="int", help="Size of shift")
    p.add_option("--subtract", help="Subtract bases from window")
    p.add_option(
        "--nomerge", default=False, action="store_true", help="Do not merge features"
    )


def check_window_options(opts):
    window = opts.window
    shift = opts.shift
    subtract = opts.subtract
    assert window % shift == 0, "--window must be divisible by --shift"
    logging.debug(
        "Line/stack-plot options: window={0} shift={1} subtract={2}".format(
            window, shift, subtract
        )
    )
    merge = not opts.nomerge

    return window, shift, subtract, merge


def get_beds(s, binned=False):
    return [x + ".bed" for x in s] if not binned else [x for x in s]


def get_nbins(clen, shift):
    nbins = clen / shift
    if clen % shift:
        nbins += 1
    return nbins


def linearray(binfile, chr, window, shift):
    mn = binfile.mapping[chr]
    m, n = zip(*mn)

    m = np.array(m, dtype="float")
    w = window / shift
    m = moving_sum(m, window=w)
    return m


def lineplot(ax, binfiles, nbins, chr, window, shift, color="br"):
    assert len(binfiles) <= 2, "A max of two line plots are supported"

    t = np.arange(nbins)
    bf = binfiles[0]
    m = linearray(bf, chr, window, shift)
    ax.plot(t, m, "{0}-".format(color[0]), lw=2)

    formatter = ticker.FuncFormatter(
        lambda x, pos: human_readable_base(int(x) * shift, pos)
    )
    ax.xaxis.set_major_formatter(formatter)
    for tl in ax.get_xticklabels():
        tl.set_color("darkslategray")

    label = bf.filename.split(".")[0]
    perw = "per {0}".format(human_size(window, precision=0))
    ax.set_ylabel(label + " " + perw, color=color[0])

    if len(binfiles) == 2:
        ax2 = ax.twinx()
        bf = binfiles[1]
        m = linearray(bf, chr, window, shift)
        ax2.plot(t, m, "{0}-".format(color[1]), lw=2)
        # Differentiate tick labels through colors
        for tl in ax.get_yticklabels():
            tl.set_color(color[0])
        for tl in ax2.get_yticklabels():
            tl.set_color(color[1])

        label = bf.filename.split(".")[0]
        ax2.set_ylabel(label + " " + perw, color=color[1])

    ax.set_xlim(0, nbins)


def composite(args):
    """
    %prog composite fastafile chr1

    Combine line plots, feature bars and alt-bars, different data types
    specified in options. Inputs must be BED-formatted. Three types of viz are
    currently supported:

    --lines: traditional line plots, useful for plotting feature freq
    --bars: show where the extent of features are
    --altbars: similar to bars, yet in two alternating tracks, e.g. scaffolds
    """
    from jcvi.graphics.chromosome import HorizontalChromosome

    p = OptionParser(composite.__doc__)
    p.add_option("--lines", help="Features to plot in lineplot")
    p.add_option("--bars", help="Features to plot in bars")
    p.add_option("--altbars", help="Features to plot in alt-bars")
    p.add_option(
        "--fatten",
        default=False,
        action="store_true",
        help="Help visualize certain narrow features",
    )
    p.add_option(
        "--mode",
        default="span",
        choices=("span", "count", "score"),
        help="Accumulate feature based on",
    )
    add_window_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="8x5")

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, chr = args
    window, shift, subtract, merge = check_window_options(opts)
    linebeds, barbeds, altbarbeds = [], [], []
    fatten = opts.fatten
    if opts.lines:
        lines = opts.lines.split(",")
        linebeds = get_beds(lines)
    if opts.bars:
        bars = opts.bars.split(",")
        barbeds = get_beds(bars)
    if opts.altbars:
        altbars = opts.altbars.split(",")
        altbarbeds = get_beds(altbars)

    linebins = get_binfiles(linebeds, fastafile, shift, mode=opts.mode, merge=merge)

    margin = 0.12
    clen = Sizes(fastafile).mapping[chr]
    nbins = get_nbins(clen, shift)

    plt.rcParams["xtick.major.size"] = 0
    plt.rcParams["ytick.major.size"] = 0

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    root.text(0.5, 0.95, chr, ha="center", color="darkslategray")

    xstart, xend = margin, 1 - margin
    xlen = xend - xstart
    ratio = xlen / clen
    # Line plots
    ax = fig.add_axes([xstart, 0.6, xlen, 0.3])
    lineplot(ax, linebins, nbins, chr, window, shift)

    # Bar plots
    yy = 0.5
    yinterval = 0.08
    xs = lambda x: xstart + ratio * x
    r = 0.01
    fattend = 0.0025
    for bb in barbeds:
        root.text(xend + 0.01, yy, bb.split(".")[0], va="center")
        HorizontalChromosome(root, xstart, xend, yy, height=0.02)
        bb = Bed(bb)
        for b in bb:
            start, end = xs(b.start), xs(b.end)
            span = end - start
            if fatten and span < fattend:
                span = fattend

            root.add_patch(
                Rectangle((start, yy - r), span, 2 * r, lw=0, fc="darkslategray")
            )
        yy -= yinterval

    # Alternative bar plots
    offset = r / 2
    for bb in altbarbeds:
        root.text(xend + 0.01, yy, bb.split(".")[0], va="center")
        bb = Bed(bb)
        for i, b in enumerate(bb):
            start, end = xs(b.start), xs(b.end)
            span = end - start
            if span < 0.0001:
                continue
            offset = -offset
            root.add_patch(
                Rectangle(
                    (start, yy + offset), end - start, 0.003, lw=0, fc="darkslategray"
                )
            )
        yy -= yinterval

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = chr + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def multilineplot(args):
    """
    %prog multilineplot fastafile chr1

    Combine multiple line plots in one vertical stack
    Inputs must be BED-formatted.

    --lines: traditional line plots, useful for plotting feature freq
    """
    p = OptionParser(multilineplot.__doc__)
    p.add_option("--lines", help="Features to plot in lineplot")
    p.add_option("--colors", help="List of colors matching number of input bed files")
    p.add_option(
        "--mode",
        default="span",
        choices=("span", "count", "score"),
        help="Accumulate feature based on",
    )
    p.add_option(
        "--binned",
        default=False,
        action="store_true",
        help="Specify whether the input is already binned; "
        + "if True, input files are considered to be binfiles",
    )
    p.add_option("--ymax", type="int", help="Set Y-axis max")
    add_window_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="8x5")

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, chr = args
    window, shift, subtract, merge = check_window_options(opts)
    linebeds = []
    colors = opts.colors
    if opts.lines:
        lines = opts.lines.split(",")
        assert len(colors) == len(lines), (
            "Number of chosen colors must match" + " number of input bed files"
        )
        linebeds = get_beds(lines, binned=opts.binned)

    linebins = get_binfiles(
        linebeds, fastafile, shift, mode=opts.mode, binned=opts.binned, merge=merge
    )

    clen = Sizes(fastafile).mapping[chr]
    nbins = get_nbins(clen, shift)

    plt.rcParams["xtick.major.size"] = 0
    plt.rcParams["ytick.major.size"] = 0
    plt.rcParams["figure.figsize"] = iopts.w, iopts.h

    fig, axarr = plt.subplots(nrows=len(lines))
    if len(linebeds) == 1:
        axarr = (axarr,)
    fig.suptitle(latex(chr), color="darkslategray")

    for i, ax in enumerate(axarr):
        lineplot(
            ax,
            [linebins[i]],
            nbins,
            chr,
            window,
            shift,
            color="{0}{1}".format(colors[i], "r"),
        )

    if opts.ymax:
        ax.set_ylim(0, opts.ymax)

    plt.subplots_adjust(hspace=0.5)

    image_name = chr + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def heatmap(args):
    """
    %prog heatmap fastafile chr1

    Combine stack plot with heatmap to show abundance of various tracks along
    given chromosome. Need to give multiple beds to --stacks and --heatmaps
    """
    p = OptionParser(heatmap.__doc__)
    p.add_option(
        "--stacks",
        default="Exons,Introns,DNA_transposons,Retrotransposons",
        help="Features to plot in stackplot",
    )
    p.add_option(
        "--heatmaps",
        default="Copia,Gypsy,hAT,Helitron,Introns,Exons",
        help="Features to plot in heatmaps",
    )
    p.add_option("--meres", default=None, help="Extra centromere / telomere features")
    add_window_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="8x5")

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, chr = args
    window, shift, subtract, merge = check_window_options(opts)

    stacks = opts.stacks.split(",")
    heatmaps = opts.heatmaps.split(",")
    stackbeds = get_beds(stacks)
    heatmapbeds = get_beds(heatmaps)
    stackbins = get_binfiles(
        stackbeds, fastafile, shift, subtract=subtract, merge=merge
    )
    heatmapbins = get_binfiles(
        heatmapbeds, fastafile, shift, subtract=subtract, merge=merge
    )

    margin = 0.06
    inner = 0.015
    clen = Sizes(fastafile).mapping[chr]

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Gauge
    ratio = draw_gauge(root, margin, clen, rightmargin=4 * margin)
    yinterval = 0.3
    xx = margin
    yy = 1 - margin
    yy -= yinterval
    xlen = clen / ratio
    cc = chr
    if "_" in chr:
        ca, cb = chr.split("_")
        cc = ca[0].upper() + cb

    root.add_patch(Rectangle((xx, yy), xlen, yinterval - inner, color=gray))
    ax = fig.add_axes([xx, yy, xlen, yinterval - inner])

    nbins = get_nbins(clen, shift)

    owindow = clen / 100
    if owindow > window:
        window = owindow / shift * shift

    stackplot(ax, stackbins, nbins, palette, chr, window, shift)
    ax.text(
        0.1,
        0.9,
        cc,
        va="top",
        zorder=100,
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", fc="w", alpha=0.5),
    )

    # Legends
    xx += xlen + 0.01
    yspace = (yinterval - inner) / (len(stackbins) + 1)
    yy = 1 - margin - yinterval
    for s, p in zip(stacks, palette):
        s = s.replace("_", " ")
        s = Registration.get(s, s)

        yy += yspace
        root.add_patch(Rectangle((xx, yy), inner, inner, color=p, lw=0))
        root.text(xx + 1.5 * inner, yy, s, size=10)

    yh = 0.05  # Heatmap height
    # Heatmaps
    xx = margin
    yy = 1 - margin - yinterval - inner
    for s, p in zip(heatmaps, heatmapbins):
        s = s.replace("_", " ")
        s = Registration.get(s, s)

        yy -= yh
        m = stackarray(p, chr, window, shift)

        Y = np.array([m, m])
        root.imshow(
            Y,
            extent=(xx, xx + xlen, yy, yy + yh - inner),
            interpolation="nearest",
            aspect="auto",
            cmap=iopts.cmap,
        )
        root.text(xx + xlen + 0.01, yy, s, size=10)

    yy -= yh

    meres = opts.meres
    if meres:
        bed = Bed(meres)
        for b in bed:
            if b.seqid != chr:
                continue
            pos = (b.start + b.end) / 2
            cpos = pos / ratio
            xx = margin + cpos
            accn = b.accn.capitalize()
            root.add_patch(CirclePolygon((xx, yy), radius=0.01, fc="m", ec="m"))
            root.text(xx + 0.014, yy, accn, va="center", color="m")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = chr + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def draw_gauge(ax, margin, maxl, rightmargin=None):
    # Draw a gauge on the top of the canvas
    rightmargin = rightmargin or margin
    ax.plot([margin, 1 - rightmargin], [1 - margin, 1 - margin], "k-", lw=2)

    best_stride = autoscale(maxl)
    nintervals = maxl * 1.0 / best_stride

    xx, yy = margin, 1 - margin
    tip = 0.005
    xinterval = (1 - margin - rightmargin) / nintervals
    l = human_size(best_stride)
    if l[-1] == "b":
        suffix = target = l[-2:]

    for i in range(0, maxl + 1, best_stride):
        l = human_size(i, precision=0, target=target)
        if l[-1] == "b":
            l, suffix = l[:-2], l[-2:]
        ax.plot([xx, xx], [yy, yy + tip], "k-", lw=2)
        ax.text(xx, yy + 2 * tip, l, ha="center", size=13)
        xx += xinterval

    xx += 4 * tip - xinterval
    ax.text(xx + tip, yy + 2 * tip, suffix)

    return best_stride / xinterval


def get_binfiles(
    inputfiles, fastafile, shift, mode="span", subtract=None, binned=False, merge=True
):
    if not binned:
        binopts = ["--binsize={0}".format(shift)]
        binopts.append("--mode={0}".format(mode))
        if subtract:
            binopts.append("--subtract={0}".format(subtract))
        if not merge:
            binopts.append("--nomerge")
        binfiles = [bins([x, fastafile] + binopts) for x in inputfiles if op.exists(x)]
    else:
        binfiles = inputfiles
    binfiles = [BinFile(x) for x in binfiles]

    return binfiles


def stackarray(binfile, chr, window, shift):
    mn = binfile.mapping[chr]
    m, n = zip(*mn)

    m = np.array(m, dtype="float")
    n = np.array(n, dtype="float")

    w = window / shift
    m = moving_sum(m, window=w)
    n = moving_sum(n, window=w)
    m /= n

    return m


def stackplot(ax, binfiles, nbins, palette, chr, window, shift):
    t = np.arange(nbins, dtype="float") + 0.5
    m = np.zeros(nbins, dtype="float")
    zorders = range(10)[::-1]
    for binfile, p, z in zip(binfiles, palette, zorders):
        s = stackarray(binfile, chr, window, shift)
        m += s
        ax.fill_between(t, m, color=p, lw=0, zorder=z)

    ax.set_xlim(0, nbins)
    ax.set_ylim(0, 1)
    ax.set_axis_off()


def stack(args):
    """
    %prog stack fastafile

    Create landscape plots that show the amounts of genic sequences, and repetitive
    sequences along the chromosomes.
    """
    p = OptionParser(stack.__doc__)
    p.add_option("--top", default=10, type="int", help="Draw the first N chromosomes")
    p.add_option(
        "--stacks",
        default="Exons,Introns,DNA_transposons,Retrotransposons",
        help="Features to plot in stackplot",
    )
    p.add_option("--switch", help="Change chr names based on two-column file")
    add_window_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="8x8")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    top = opts.top
    window, shift, subtract, merge = check_window_options(opts)
    switch = opts.switch
    if switch:
        switch = DictFile(opts.switch)

    stacks = opts.stacks.split(",")
    bedfiles = get_beds(stacks)
    binfiles = get_binfiles(bedfiles, fastafile, shift, subtract=subtract, merge=merge)

    sizes = Sizes(fastafile)
    s = list(sizes.iter_sizes())[:top]
    maxl = max(x[1] for x in s)
    margin = 0.08
    inner = 0.02  # y distance between tracks

    pf = fastafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Gauge
    ratio = draw_gauge(root, margin, maxl)

    # Per chromosome
    yinterval = (1 - 2 * margin) / (top + 1)
    xx = margin
    yy = 1 - margin
    for chr, clen in s:
        yy -= yinterval
        xlen = clen / ratio
        cc = chr
        if "_" in chr:
            ca, cb = chr.split("_")
            cc = ca[0].upper() + cb

        if switch and cc in switch:
            cc = "\n".join((cc, "({0})".format(switch[cc])))

        root.add_patch(Rectangle((xx, yy), xlen, yinterval - inner, color=gray))
        ax = fig.add_axes([xx, yy, xlen, yinterval - inner])

        nbins = clen / shift
        if clen % shift:
            nbins += 1

        stackplot(ax, binfiles, nbins, palette, chr, window, shift)
        root.text(
            xx - 0.04, yy + 0.5 * (yinterval - inner), cc, ha="center", va="center"
        )

        ax.set_xlim(0, nbins)
        ax.set_ylim(0, 1)
        ax.set_axis_off()

    # Legends
    yy -= yinterval
    xx = margin
    for b, p in zip(bedfiles, palette):
        b = b.rsplit(".", 1)[0].replace("_", " ")
        b = Registration.get(b, b)

        root.add_patch(Rectangle((xx, yy), inner, inner, color=p, lw=0))
        xx += 2 * inner
        root.text(xx, yy, b, size=13)
        xx += len(b) * 0.012 + inner

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main()
