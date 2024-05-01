#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Create chromosome landscape plots that are similar to the ones used in soybean
and sorghum paper.
"""


import os.path as op
import sys

from collections import Counter, OrderedDict, defaultdict
from typing import List, Optional

import numpy as np

from ..algorithms.matrix import moving_sum
from ..apps.base import ActionDispatcher, OptionParser, logger
from ..formats.base import BaseFile, DictFile, LineFile, must_open
from ..formats.bed import Bed, bins, get_nbins
from ..formats.sizes import Sizes
from ..utils.cbook import autoscale, human_size

from .base import (
    CirclePolygon,
    Colormap,
    Extent,
    Rectangle,
    adjust_extent,
    adjust_spines,
    human_readable_base,
    latex,
    markup,
    normalize_axes,
    plt,
    savefig,
    set_human_axis,
    ticker,
)
from .chromosome import HorizontalChromosome

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


class BinLine:
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

        fp = open(filename, encoding="utf-8")
        for row in fp:
            b = BinLine(row)
            self.append(b)
            chr, len, binlen = b.chr, b.len, b.binlen
            self.mapping[chr].append((len, binlen))
        fp.close()


class ChrInfoLine:
    def __init__(self, row, delimiter=","):
        args = [x.strip() for x in row.split(delimiter)]
        self.name = args[0]
        self.color = args[1]
        if len(args) > 2:
            self.new_name = args[2]
        else:
            self.new_name = self.name


class ChrInfoFile(BaseFile, OrderedDict):
    def __init__(self, filename, delimiter=","):
        super(ChrInfoFile, self).__init__(filename)
        with open(filename, encoding="utf-8") as fp:
            for row in fp:
                if row[0] == "#":
                    continue
                line = ChrInfoLine(row, delimiter=delimiter)
                self[line.name] = line


class TitleInfoLine:
    def __init__(self, row, delimiter=","):
        args = [x.strip() for x in row.split(delimiter)]
        self.name = args[0]
        self.title = args[1]
        self.subtitle = None
        if len(args) > 2:
            self.subtitle = args[2]


class TitleInfoFile(BaseFile, OrderedDict):
    def __init__(self, filename, delimiter=","):
        super(TitleInfoFile, self).__init__(filename)
        with open(filename, encoding="utf-8") as fp:
            for row in fp:
                if row[0] == "#":
                    continue
                line = TitleInfoLine(row, delimiter=delimiter)
                self[line.name] = line


def main():

    actions = (
        ("composite", "combine line plots, feature bars and alt-bars"),
        ("depth", "show per chromosome depth plot across genome"),
        ("heatmap", "similar to stack but adding heatmap"),
        ("mosdepth", "plot depth vs. coverage per chromosome"),
        ("multilineplot", "combine multiple line plots in one vertical stack"),
        ("stack", "create landscape plot with genic/te composition"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def parse_distfile(filename):
    """Parse mosdepth dist.txt file. The file has contents like:

    #chr    start   end     depth   (header added here for clarity)
    chr01A  0       50000   31.00
    chr01A  50000   100000  36.00
    chr01A  100000  150000  280.00
    chr01A  150000  200000  190.00

    Args:
        filename (str): Path to the file.
    """
    dists = defaultdict(Counter)
    with must_open(filename) as fp:
        for row in fp:
            chromosome, _, _, depth = row.split()
            depth = int(float(depth))
            dists[chromosome][depth] += 1
    logger.debug("Loaded %d seqids", len(dists))
    return dists


def parse_groupsfile(filename):
    """Parse groupsfile, which contains the tracks to be plotted
    in the vertically stacked mosdepth plot.

    chr01A,chr01B g,m
    chr02A,chr02B g,m
    chr03A,chr03B g,m

    Args:
        filename (str): Path to the groups file.
    """
    groups = []
    with open(filename, encoding="utf-8") as fp:
        for row in fp:
            chrs, colors = row.split()
            groups.append((chrs.split(","), colors.split(",")))
    logger.debug("Loaded %d groups", len(groups))
    return groups


def cumarray_to_array(ar):
    """Convert cumulative array to normal array.

    Args:
        ar (List): List of numbers
    """
    ans = []
    for i, x in enumerate(ar):
        ans.append(x if i == 0 else (ar[i] - ar[i - 1]))
    return ans


def mosdepth(args):
    """
    %prog mosdepth mosdepth.global.dist.txt groups

    Plot depth vs. coverage per chromosome. Inspired by mosdepth plot. See also:
    https://github.com/brentp/mosdepth
    """
    import seaborn as sns

    sns.set_style("darkgrid")

    p = OptionParser(mosdepth.__doc__)
    p.add_argument("--maxdepth", default=100, type=int, help="Maximum depth to plot")
    p.add_argument(
        "--logscale", default=False, action="store_true", help="Use log-scale on depth"
    )
    opts, args, iopts = p.set_image_options(args, style="dark", figsize="6x8")

    if len(args) != 2:
        sys.exit(p.print_help())

    # Read in datasets
    distfile, groupsfile = args
    dists = parse_distfile(distfile)
    groups = parse_groupsfile(groupsfile)
    logscale = opts.logscale

    # Construct a composite figure with N tracks indicated in the groups
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))

    rows = len(groups)
    ypad = 0.05
    yinterval = (1 - 2 * ypad) / (rows + 1)
    yy = 1 - ypad

    for group_idx, (chrs, colors) in enumerate(groups):
        yy -= yinterval
        ax = fig.add_axes((0.15, yy, 0.7, yinterval * 0.85))
        for c, color in zip(chrs, colors):
            cdata = dists[c].items()
            logger.debug("Importing %d records for %s", len(cdata), c)
            cx, cy = zip(*sorted(cdata))
            ax.plot(cx, cy, "-", color=color)
        if logscale:
            ax.set_xscale("log", basex=2)
        ax.set_xlim(1 if logscale else 0, opts.maxdepth)
        ax.get_yaxis().set_visible(False)
        if group_idx != rows - 1:
            ax.get_xaxis().set_visible(False)

        # Add legend to the right of the canvas
        label_pad = 0.02
        label_yy = yy + yinterval
        for c, color in zip(chrs, colors):
            label_yy -= label_pad
            root.text(0.92, label_yy, c, color=color, ha="center", va="center")

    root.text(
        0.1,
        0.5,
        "Proportion of bases at coverage",
        rotation=90,
        color="darkslategray",
        ha="center",
        va="center",
    )
    root.text(0.5, 0.05, "Coverage", color="darkslategray", ha="center", va="center")
    normalize_axes(root)
    adjust_spines(ax, ["bottom"], outward=True)

    pf = "mosdepth"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def draw_depth(
    root,
    ax,
    bed: Bed,
    chrinfo: dict = {},
    defaultcolor: str = "k",
    sepcolor: str = "w",
    maxdepth: int = 100,
    logscale: bool = False,
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
):
    """Draw depth plot on the given axes, using data from bed

    Args:
        root (matplotlib.Axes): Canvas axes
        ax (matplotlib.Axes): Axes to plot data on
        bed (Bed): Bed data from mosdepth
        chrinfo (ChrInfoFile): seqid => color, new name
        defaultcolor (str): matplotlib-compatible color for data points
        sepcolor (str): matplotlib-compatible color for chromosome breaks
        maxdepth (int): Upper limit of the y-axis (depth)
        title (str): Title of the figure, to the right of the axis
        subtitle (str): Subtitle of the figure, just below title
    """
    if chrinfo is None:
        chrinfo = {}
    sizes = bed.max_bp_in_chr
    seqids = chrinfo.keys() if chrinfo else sizes.keys()
    starts = {}
    ends = {}
    label_positions = []
    start = 0
    for seqid in seqids:
        if seqid not in sizes:
            continue
        starts[seqid] = start
        end = start + sizes[seqid]
        ends[seqid] = end
        label_positions.append((seqid, (start + end) / 2))
        start = end
    xsize = end

    # Extract plotting data
    data = []
    data_by_seqid = defaultdict(list)
    for b in bed:
        seqid = b.seqid
        if seqid not in starts:
            continue
        # chr01A  2000000 3000000 113.00
        x = starts[seqid] + (b.start + b.end) / 2
        y = float(b.accn)
        c = chrinfo[seqid].color if seqid in chrinfo else "k"
        data.append((x, y, c))
        data_by_seqid[seqid].append(y)

    x, y, c = zip(*data)
    ax.scatter(
        x,
        y,
        c=c,
        edgecolors="none",
        s=8,
        lw=0,
    )
    logger.debug("Obtained %d data points with depth data", len(data))

    # Per seqid median
    medians = {}
    for seqid, values in data_by_seqid.items():
        c = chrinfo[seqid].color if seqid in chrinfo else defaultcolor
        seqid_start = starts[seqid]
        seqid_end = ends[seqid]
        seqid_median = np.median(values)
        medians[seqid] = seqid_median
        ax.plot(
            (seqid_start, seqid_end),
            (seqid_median, seqid_median),
            "-",
            lw=4,
            color=c,
            alpha=0.5,
        )

    # Vertical lines for all the breaks
    for pos in starts.values():
        ax.plot((pos, pos), (0, maxdepth), "-", lw=1, color=sepcolor)

    # Beautify the numeric axis
    for tick in ax.get_xticklines() + ax.get_yticklines():
        tick.set_visible(False)

    median_depth_y = 0.88
    chr_label_y = 0.08
    for seqid, position in label_positions:
        xpos = 0.1 + position * 0.8 / xsize
        c = chrinfo[seqid].color if seqid in chrinfo else defaultcolor
        newseqid = chrinfo[seqid].new_name if seqid in chrinfo else seqid
        root.text(
            xpos, chr_label_y, newseqid, color=c, ha="center", va="center", rotation=20
        )
        seqid_median = medians[seqid]
        root.text(
            xpos,
            median_depth_y,
            str(int(seqid_median)),
            color=c,
            ha="center",
            va="center",
        )

    # Add an arrow to the right of the plot, indicating these are median depths
    root.text(
        0.91,
        0.88,
        r"$\leftarrow$median",
        color="lightslategray",
        va="center",
    )

    if title:
        root.text(
            0.95,
            0.5,
            markup(title),
            color="darkslategray",
            ha="center",
            va="center",
            size=15,
        )
    if subtitle:
        root.text(
            0.95,
            0.375,
            markup(subtitle),
            color="darkslategray",
            ha="center",
            va="center",
            size=15,
        )

    ax.set_xticks([])
    ax.set_xlim(0, xsize)
    if logscale:
        ax.set_yscale("log", basey=2)
    ax.set_ylim(1 if logscale else 0, maxdepth)
    ax.set_ylabel("Depth")

    set_human_axis(ax)
    plt.setp(ax.get_xticklabels() + ax.get_yticklabels(), color="gray", size=10)
    normalize_axes(root)


def draw_multi_depth(
    root,
    panel_roots,
    panel_axes,
    bedfiles: List[str],
    chrinfo_file: str,
    titleinfo_file: str,
    maxdepth: int,
    logscale: bool,
):
    """
    Draw multiple depth plots on the same canvas.
    """
    chrinfo = ChrInfoFile(chrinfo_file) if chrinfo_file else {}
    titleinfo = TitleInfoFile(titleinfo_file) if titleinfo_file else {}
    npanels = len(bedfiles)
    yinterval = 1.0 / npanels
    ypos = 1 - yinterval
    for bedfile, panel_root, panel_ax in zip(bedfiles, panel_roots, panel_axes):
        pf = op.basename(bedfile).split(".", 1)[0]
        bed = Bed(bedfile)

        if ypos > 0.001:
            root.plot((0.02, 0.98), (ypos, ypos), "-", lw=2, color="lightslategray")

        title = titleinfo.get(bedfile, pf.split("_", 1)[0])
        subtitle = None
        if isinstance(title, TitleInfoLine):
            subtitle = title.subtitle
            title = title.title

        draw_depth(
            panel_root,
            panel_ax,
            bed,
            chrinfo=chrinfo,
            maxdepth=maxdepth,
            logscale=logscale,
            title=title,
            subtitle=subtitle,
        )
        ypos -= yinterval

    normalize_axes(root)


def depth(args):
    """
    %prog depth *.regions.bed.gz

    Plot the mosdepth regions BED file. We recommend to generate this BED file
    by (please adjust the --by parameter to your required resolution):

    $ mosdepth --no-per-base --use-median --fast-mode --by 1000000 sample.wgs
    sample.bam

    Use --chrinfo to specify a colormap between seqid, desired color, and
    optionally a new name. For example:

    chr01A, #c51b7d, 1A
    chr01B, #4d9221, 1B
    ...

    Only seqids that are in the colormap will be plotted, in the order that's
    given in the file. When --colormap is not set, every seqid will be drawn in
    black.

    Can take multiple BED files as input and then plot all of them in a
    composite figure.
    """
    p = OptionParser(depth.__doc__)
    p.add_argument(
        "--chrinfo", help="Comma-separated mappings between seqid, color, new_name"
    )
    p.add_argument(
        "--titleinfo",
        help="Comma-separated titles mappings between filename, title",
    )
    p.add_argument("--maxdepth", default=100, type=int, help="Maximum depth to show")
    p.add_argument(
        "--logscale", default=False, action="store_true", help="Use log-scale on depth"
    )
    opts, args, iopts = p.set_image_options(args, style="dark", figsize="14x4")

    if len(args) < 1:
        sys.exit(not p.print_help())

    bedfiles = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))

    npanels = len(bedfiles)
    yinterval = 1.0 / npanels
    ypos = 1 - yinterval
    panel_roots, panel_axes = [], []
    for _ in range(npanels):
        panel_root = root if npanels == 1 else fig.add_axes((0, ypos, 1, yinterval))
        panel_ax = fig.add_axes((0.1, ypos + 0.2 * yinterval, 0.8, 0.65 * yinterval))
        panel_roots.append(panel_root)
        panel_axes.append(panel_ax)
        ypos -= yinterval

    draw_multi_depth(
        root,
        panel_roots,
        panel_axes,
        bedfiles,
        opts.chrinfo,
        opts.titleinfo,
        opts.maxdepth,
        opts.logscale,
    )

    if npanels > 1:
        pf = op.commonprefix(bedfiles)
    pf = pf or "depth"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def add_window_options(p):
    """
    Add options for window plotting.
    """
    p.add_argument("--window", default=500000, type=int, help="Size of window")
    p.add_argument("--shift", default=100000, type=int, help="Size of shift")
    p.add_argument("--subtract", help="Subtract bases from window")
    p.add_argument(
        "--nomerge", default=False, action="store_true", help="Do not merge features"
    )


def check_window_options(opts):
    """
    Check the window options, and return the values.
    """
    window = opts.window
    shift = opts.shift
    subtract = opts.subtract
    assert window % shift == 0, "--window must be divisible by --shift"
    logger.debug(
        "Line/stack-plot options: window=%d shift=%d subtract=%s",
        window,
        shift,
        subtract,
    )
    merge = not opts.nomerge

    return window, shift, subtract, merge


def get_beds(s: List[str], binned: bool = False) -> List[str]:
    """
    Get the bed files for each feature, and return them as a list.
    """
    return [x + ".bed" for x in s] if not binned else [x for x in s]


def linearray(binfile, chr, window, shift):
    mn = binfile.mapping[chr]
    m, _ = zip(*mn)

    m = np.array(m, dtype=float)
    w = window // shift
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
    p = OptionParser(composite.__doc__)
    p.add_argument("--lines", help="Features to plot in lineplot")
    p.add_argument("--bars", help="Features to plot in bars")
    p.add_argument("--altbars", help="Features to plot in alt-bars")
    p.add_argument(
        "--fatten",
        default=False,
        action="store_true",
        help="Help visualize certain narrow features",
    )
    p.add_argument(
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
    window, shift, _, merge = check_window_options(opts)
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
    nbins, _ = get_nbins(clen, shift)

    plt.rcParams["xtick.major.size"] = 0
    plt.rcParams["ytick.major.size"] = 0

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))

    root.text(0.5, 0.95, chr, ha="center", color="darkslategray")

    xstart, xend = margin, 1 - margin
    xlen = xend - xstart
    ratio = xlen / clen
    # Line plots
    ax = fig.add_axes((xstart, 0.6, xlen, 0.3))
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
        for b in bb:
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
    p.add_argument("--lines", help="Features to plot in lineplot")
    p.add_argument("--colors", help="List of colors matching number of input bed files")
    p.add_argument(
        "--mode",
        default="span",
        choices=("span", "count", "score"),
        help="Accumulate feature based on",
    )
    p.add_argument(
        "--binned",
        default=False,
        action="store_true",
        help="Specify whether the input is already binned; "
        + "if True, input files are considered to be binfiles",
    )
    p.add_argument("--ymax", type=int, help="Set Y-axis max")
    add_window_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="8x5")

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, chr = args
    window, shift, _, merge = check_window_options(opts)
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
    nbins, _ = get_nbins(clen, shift)

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


def draw_heatmaps(
    fig,
    root,
    root_extent: Extent,
    fastafile: str,
    chr: str,
    stacks: List[str],
    heatmaps: List[str],
    window: int,
    shift: int,
    cmap: Colormap,
    subtract: Optional[int] = None,
    merge: bool = False,
    meres: Optional[str] = None,
):
    """
    Draw heatmap for the given chromosome.
    """
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
    extent = (xx, yy, xlen, yinterval - inner)
    adjusted = adjust_extent(extent, root_extent)
    ax = fig.add_axes(adjusted)

    nbins, _ = get_nbins(clen, shift)

    owindow = clen / 100
    if owindow > window:
        window = owindow // shift * shift

    stackplot(ax, stackbins, nbins, palette, chr, window, shift)
    ax.text(
        0.05,
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
            cmap=cmap,
        )
        root.text(xx + xlen + 0.01, yy, s, size=10)

    yy -= yh

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

    normalize_axes(root)


def heatmap(args):
    """
    %prog heatmap fastafile chr1

    Combine stack plot with heatmap to show abundance of various tracks along
    given chromosome. Need to give multiple beds to --stacks and --heatmaps
    """
    p = OptionParser(heatmap.__doc__)
    p.add_argument(
        "--stacks",
        default="Exons,Introns,DNA_transposons,Retrotransposons",
        help="Features to plot in stackplot",
    )
    p.add_argument(
        "--heatmaps",
        default="Copia,Gypsy,hAT,Helitron,Introns,Exons",
        help="Features to plot in heatmaps",
    )
    p.add_argument("--meres", default=None, help="Extra centromere / telomere features")
    add_window_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="8x5")

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, chr = args
    window, shift, subtract, merge = check_window_options(opts)

    stacks = opts.stacks.split(",")
    heatmaps = opts.heatmaps.split(",")

    fig = plt.figure(1, (iopts.w, iopts.h))
    root_extent = (0, 0, 1, 1)
    root = fig.add_axes(root_extent)

    draw_heatmaps(
        fig,
        root,
        root_extent,
        fastafile,
        chr,
        stacks,
        heatmaps,
        window,
        shift,
        iopts.cmap,
        subtract,
        merge,
        meres=opts.meres,
    )

    image_name = chr + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def draw_gauge(ax, margin: float, maxl: int, rightmargin: Optional[float] = None):
    """
    Draw a gauge on the top of the canvas, showing the scale of the chromosome.
    """
    rightmargin = rightmargin or margin
    ax.plot([margin, 1 - rightmargin], [1 - margin, 1 - margin], "k-", lw=2)

    best_stride = autoscale(maxl)
    nintervals = maxl / best_stride

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
    inputfiles: List[str],
    fastafile: str,
    shift: int,
    mode: str = "span",
    subtract: Optional[int] = None,
    binned: bool = False,
    merge: bool = True,
):
    """
    Get binfiles from input files. If not binned, then bin them first.
    """
    if not binned:
        binopts = [f"--binsize={shift}"]
        binopts.append(f"--mode={mode}")
        if subtract:
            binopts.append(f"--subtract={subtract}")
        if not merge:
            binopts.append("--nomerge")
        binfiles = [bins([x, fastafile] + binopts) for x in inputfiles if op.exists(x)]
    else:
        binfiles = inputfiles
    binfiles = [BinFile(x) for x in binfiles]

    return binfiles


def stackarray(binfile: BinFile, chr: str, window: int, shift: int):
    """
    Get stack array from binfile for the given chr.
    """
    mn = binfile.mapping[chr]
    m, n = zip(*mn)

    m = np.array(m, dtype=float)
    n = np.array(n, dtype=float)

    w = window // shift
    m = moving_sum(m, window=w)
    n = moving_sum(n, window=w)
    m /= n

    return m


def stackplot(
    ax,
    binfiles: List[BinFile],
    nbins: int,
    palette: List[str],
    chr: str,
    window: int,
    shift: int,
):
    """
    Plot stackplot on the given axes, using data from binfiles.
    """
    t = np.arange(nbins, dtype=float) + 0.5
    m = np.zeros(nbins, dtype=float)
    zorders = range(10)[::-1]
    for binfile, p, z in zip(binfiles, palette, zorders):
        s = stackarray(binfile, chr, window, shift)
        m += s
        ax.fill_between(t, m, color=p, lw=0, zorder=z)

    ax.set_xlim(0, nbins)
    ax.set_ylim(0, 1)
    ax.set_axis_off()


def draw_stacks(
    fig,
    root,
    root_extent: Extent,
    stacks: List[str],
    fastafile: str,
    window: int,
    shift: int,
    top: int,
    merge: bool = True,
    subtract: Optional[int] = None,
    switch: Optional[DictFile] = None,
):
    """
    Draw stack plot.
    """
    bedfiles = get_beds(stacks)
    binfiles = get_binfiles(bedfiles, fastafile, shift, subtract=subtract, merge=merge)

    sizes = Sizes(fastafile)
    s = list(sizes.iter_sizes())[:top]
    maxl = max(x[1] for x in s)
    margin = 0.08
    inner = 0.02  # y distance between tracks

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
            cc = "\n".join((cc, f"({switch[cc]})"))

        extent = (xx, yy, xlen, yinterval - inner)
        adjusted = adjust_extent(extent, root_extent)
        root.add_patch(Rectangle((xx, yy), xlen, yinterval - inner, color=gray))
        ax = fig.add_axes(adjusted)

        nbins, _ = get_nbins(clen, shift)

        stackplot(ax, binfiles, nbins, palette, chr, window, shift)
        root.text(
            xx - 0.04, yy + 0.5 * (yinterval - inner), cc, ha="center", va="center"
        )

    # Legends
    yy -= yinterval
    xx = margin
    for b, p in zip(bedfiles, palette):
        b = b.rsplit(".", 1)[0].replace("_", " ")
        b = Registration.get(b, b)

        root.add_patch(Rectangle((xx, yy), inner, inner, color=p, lw=0))
        xx += 2 * inner
        root.text(xx, yy, b, size=13)
        xx += len(b) * 0.015 + inner

    normalize_axes(root)


def stack(args):
    """
    %prog stack fastafile

    Create landscape plots that show the amounts of genic sequences, and repetitive
    sequences along the chromosomes.
    """
    p = OptionParser(stack.__doc__)
    p.add_argument("--top", default=10, type=int, help="Draw the first N chromosomes")
    p.add_argument(
        "--stacks",
        default="Exons,Introns,DNA_transposons,Retrotransposons",
        help="Features to plot in stackplot",
    )
    p.add_argument("--switch", help="Change chr names based on two-column file")
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

    fig = plt.figure(1, (iopts.w, iopts.h))
    root_extent = (0, 0, 1, 1)
    root = fig.add_axes(root_extent)

    draw_stacks(
        fig,
        root,
        root_extent,
        stacks,
        fastafile,
        window,
        shift,
        top,
        merge,
        subtract,
        switch,
    )

    pf = fastafile.rsplit(".", 1)[0]
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main()
