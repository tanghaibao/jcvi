#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import copy
import os.path as op
from os import remove

import sys

from functools import partial
from typing import Optional, List, Tuple, Union

import numpy as np
import matplotlib as mpl
import seaborn as sns

mpl.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from brewer2mpl import get_map
from matplotlib import cm, rc, rcParams
from matplotlib.colors import Colormap
from matplotlib.patches import (
    Rectangle,
    Polygon,
    CirclePolygon,
    Ellipse,
    PathPatch,
    FancyArrow,
    FancyArrowPatch,
    FancyBboxPatch,
)

from ..apps.base import datadir, glob, logger, sample_N, which
from ..formats.base import LineFile
from ..utils.cbook import human_size

Extent = Tuple[float, float, float, float]

CHARS = {
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
}

GRAPHIC_FORMATS = (
    "emf",
    "eps",
    "pdf",
    "png",
    "ps",
    "raw",
    "rgba",
    "svg",
    "svgz",
)


def is_tex_available() -> bool:
    """Check if latex command is available"""
    return bool(which("latex")) and bool(which("lp"))


class ImageOptions(object):
    def __init__(self, opts):
        self.w, self.h = [int(x) for x in opts.figsize.split("x")]
        self.dpi = opts.dpi
        self.format = opts.format
        self.cmap = mpl.colormaps[opts.cmap]
        self.seed = opts.seed
        self.usetex = is_tex_available() and not opts.notex
        self.opts = opts

    def __str__(self):
        return "({0}px x {1}px)".format(self.dpi * self.w, self.dpi * self.h)

    @property
    def diverge(self):
        colors = get_map(self.opts.diverge, "diverging", 5).mpl_colors
        return colors[0], colors[-1]


class TextHandler(object):
    def __init__(self, fig, usetex: bool = True):
        self.heights = []
        try:
            self.build_height_array(fig, usetex=usetex)
        except ValueError as e:
            logger.debug(
                "Failed to init heights (error: %s). Variable label sizes skipped.", e
            )

    @classmethod
    def get_text_width_height(cls, fig, txt="chr01", size=12, usetex: bool = True):
        tp = mpl.textpath.TextPath((0, 0), txt, size=size, usetex=usetex)
        bb = tp.get_extents()
        xmin, ymin = fig.transFigure.inverted().transform((bb.xmin, bb.ymin))
        xmax, ymax = fig.transFigure.inverted().transform((bb.xmax, bb.ymax))
        return xmax - xmin, ymax - ymin

    def build_height_array(self, fig, start=1, stop=36, usetex: bool = True):
        for i in range(start, stop + 1):
            w, h = TextHandler.get_text_width_height(fig, size=i, usetex=usetex)
            self.heights.append((h, i))

    def select_fontsize(self, height, minsize=1, maxsize=12):
        if not self.heights:
            return maxsize if height > 0.01 else minsize

        from bisect import bisect_left

        i = bisect_left(self.heights, (height,))
        size = self.heights[i - 1][1] if i else minsize
        size = min(size, maxsize)
        return size


class AbstractLayout(LineFile):
    """
    Simple csv layout file for complex plotting settings. Typically, each line
    represents a subplot, a track or a panel.
    """

    def __init__(self, filename):
        super(AbstractLayout, self).__init__(filename)

    def assign_array(self, attrib, array):
        assert len(array) == len(self)
        for x, c in zip(self, array):
            if not getattr(x, attrib):
                setattr(x, attrib, c)

    def assign_colors(self, seed: Optional[int] = None):
        number = len(self)
        palette = set2_n if number <= 8 else set3_n
        # Restrict palette numbers between [3, 12]
        palette_number = max(3, min(number, 12))
        colorset = palette(palette_number)
        colorset = sample_N(colorset, number, seed=seed)
        self.assign_array("color", colorset)

    def assign_markers(self, seed: Optional[int] = None):
        markerset = sample_N(mpl.lines.Line2D.filled_markers, len(self), seed=seed)
        self.assign_array("marker", markerset)

    def __str__(self):
        return "\n".join(str(x) for x in self)


def adjust_extent(extent: Extent, root_extent: Extent) -> Extent:
    """
    Adjust the extent of the root axes.
    """
    rx, ry, rw, rh = root_extent
    ex, ey, ew, eh = extent
    return rx + ex * rw, ry + ey * rh, ew * rw, eh * rh


def linear_blend(from_color, to_color, fraction=0.5):
    """Interpolate a new color between two colors.

    https://github.com/PimpTrizkit/PJs/wiki/12.-Shade,-Blend-and-Convert-a-Web-Color-(pSBC.js)

    Args:
        from_color (matplotlib color): starting color
        to_color (matplotlib color): ending color
        fraction (float, optional): Range is 0 (closer to starting color) to 1
        (closer to ending color). Defaults to 0.5.
    """
    from matplotlib.colors import to_rgb

    def lerp(v0, v1, t):
        # Precise method, which guarantees v = v1 when t = 1
        return (1 - t) * v0 + t * v1

    r1, g1, b1 = to_rgb(from_color)
    r2, g2, b2 = to_rgb(to_color)
    return lerp(r1, r2, fraction), lerp(g1, g2, fraction), lerp(b1, b2, fraction)


def linear_shade(from_color, fraction=0.5):
    """Interpolate a lighter or darker color.

    https://github.com/PimpTrizkit/PJs/wiki/12.-Shade,-Blend-and-Convert-a-Web-Color-(pSBC.js)

    Args:
        from_color (matplotlib color): starting color
        fraction (float, optional): Range is -1 (darker) to 1 (lighter). Defaults to 0.5.
    """
    assert -1 <= fraction <= 1, "Fraction must be between -1 and 1"
    if fraction < 0:
        return linear_blend("k", from_color, 1 + fraction)
    return linear_blend(from_color, "w", fraction)


def load_image(filename: str) -> np.ndarray:
    """
    Load an image file and return as numpy array.
    """
    img = plt.imread(filename)
    if len(img.shape) == 2:  # Gray-scale image, convert to RGB
        # http://www.socouldanyone.com/2013/03/converting-grayscale-to-rgb-with-numpy.html
        h, w = img.shape
        ret = np.empty((h, w, 3), dtype=np.uint8)
        ret[:, :, 2] = ret[:, :, 1] = ret[:, :, 0] = img
        img = ret
    else:
        h, w, _ = img.shape
    logger.debug("Image `%s` loaded (%dpx x %dpx).", filename, w, h)
    return img


def latex(s):
    """Latex doesn't work well with certain characters, like '_', in plain text.
    These characters would be interpreted as control characters, so we sanitize
    these strings.

    Args:
        s (str): Input string

    Returns:
        str: Output string sanitized
    """
    return "".join([CHARS.get(char, char) for char in s])


def shorten(s, maxchar=20, mid="..."):
    if len(s) <= maxchar or len(mid) >= maxchar:
        return s
    pad = (maxchar - len(mid)) // 2
    right_pad = maxchar - len(mid) - pad
    return s[:pad] + mid + s[-right_pad:]


def set1_n(number=9):
    return get_map("Set1", "qualitative", number).hex_colors


def set2_n(number=8):
    # Get Set2 from ColorBrewer, a set of colors deemed colorblind-safe and
    # pleasant to look at by Drs. Cynthia Brewer and Mark Harrower of Pennsylvania
    # State University. These colors look lovely together, and are less
    # saturated than those colors in Set1.
    return get_map("Set2", "qualitative", number).hex_colors


def set3_n(number=12):
    return get_map("Set3", "qualitative", number).hex_colors


def paired_n(number=12):
    """See also: https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12"""
    return get_map("Paired", "qualitative", number).hex_colors


set1, set2, set3, paired = set1_n(), set2_n(), set3_n(), paired_n()


def prettyplot():
    reds = copy.copy(mpl.cm.Reds)
    reds.set_bad("white")
    reds.set_under("white")

    blues_r = copy.copy(mpl.cm.Blues_r)
    blues_r.set_bad("white")
    blues_r.set_under("white")

    # Need to 'reverse' red to blue so that blue=cold=small numbers,
    # and red=hot=large numbers with '_r' suffix
    blue_red = get_map("RdBu", "diverging", 11, reverse=True).mpl_colormap
    green_purple = get_map("PRGn", "diverging", 11).mpl_colormap
    red_purple = get_map("RdPu", "sequential", 9).mpl_colormap

    return blues_r, reds, blue_red, green_purple, red_purple


blues_r, reds, blue_red, green_purple, red_purple = prettyplot()


def normalize_axes(*axes):
    """
    Normalize the axes to have the same scale.
    """
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_axis_off()


def panel_labels(ax, labels, size: int = 16):
    """
    Add panel labels (A, B, ...) to a figure.
    """
    for xx, yy, panel_label in labels:
        if rcParams["text.usetex"]:
            panel_label = r"$\textbf{{{0}}}$".format(panel_label)
        ax.text(xx, yy, panel_label, size=size, ha="center", va="center")


def update_figname(figname: str, format: str) -> str:
    """Update the name of a figure to include the format.

    Args:
        figname (str): Path to the figure
        format (str): Figure format, must be one of GRAPHIC_FORMATS

    Returns:
        str: New file path
    """
    _, ext = op.splitext(figname)
    if ext.strip(".") in GRAPHIC_FORMATS:  # User suffix has precedence
        return figname
    # When the user has not supplied a format in the filename, use the requested format
    assert format in GRAPHIC_FORMATS, "Invalid format"
    return figname + "." + format


def savefig(figname, dpi=150, iopts=None, cleanup=True):
    try:
        format = figname.rsplit(".", 1)[-1].lower()
    except:
        format = "pdf"
    try:
        logger.debug("Matplotlib backend is: %s", mpl.get_backend())
        plt.savefig(figname, dpi=dpi, format=format)
    except Exception as e:
        logger.error("savefig failed with message:\n%s", e)
        logger.info("Try running again with --notex option to disable latex.")
        if op.exists(figname):
            if op.getsize(figname) < 1000:
                logger.debug("Cleaning up empty file: %s", figname)
                remove(figname)
        sys.exit(1)

    msg = f"Figure saved to `{figname}`"
    if iopts:
        msg += f" {iopts}"
    logger.debug(msg)

    if cleanup:
        plt.rcdefaults()


# human readable size (Kb, Mb, Gb)
def human_readable(x: Union[str, int], _, base=False):
    x = str(int(x))
    if x.endswith("000000000"):
        x = x[:-9] + "G"
    elif x.endswith("000000"):
        x = x[:-6] + "M"
    elif x.endswith("000"):
        x = x[:-3] + "K"
    if base and x[-1] in "MKG":
        x += "b"
    return x


human_readable_base = partial(human_readable, base=True)
human_formatter = ticker.FuncFormatter(human_readable)
human_base_formatter = ticker.FuncFormatter(human_readable_base)
mb_formatter = ticker.FuncFormatter(lambda x, pos: "{0}M".format(int(x / 1000000)))
mb_float_formatter = ticker.FuncFormatter(
    lambda x, pos: "{0:.1f}M".format(x / 1000000.0)
)
kb_formatter = ticker.FuncFormatter(lambda x, pos: "{0}K".format(int(x / 1000)))


def set_human_axis(ax, formatter=human_formatter):
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)


set_human_base_axis = partial(set_human_axis, formatter=human_base_formatter)


def set_helvetica_axis(ax):
    xtick_locs = ax.get_xticks().tolist()
    ytick_locs = ax.get_yticks().tolist()
    # If we dont do the following, we have
    # UserWarning: FixedFormatter should only be used together with FixedLocator
    ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(xtick_locs))
    ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(ytick_locs))
    ax.set_xticklabels([int(x) for x in xtick_locs], family="Helvetica")
    ax.set_yticklabels([int(x) for x in ytick_locs], family="Helvetica")


available_fonts = [op.basename(x) for x in glob(datadir + "/*.ttf")]


def fontprop(ax, name, size=12):
    assert name in available_fonts, "Font must be one of {0}.".format(available_fonts)

    import matplotlib.font_manager as fm

    fname = op.join(datadir, name)
    prop = fm.FontProperties(fname=fname, size=size)

    logger.debug("Set font to `%s` (`%s`)", name, prop.get_file())
    for text in ax.texts:
        text.set_fontproperties(prop)

    return prop


def markup(s):
    if "$" in s:
        return s
    import re

    s = latex(s)
    s = re.sub(r"\*(.*)\*", r"\\textit{\1}", s)
    return s


def append_percentage(s):
    # The percent symbol needs escaping in latex
    if rcParams["text.usetex"]:
        return s + r"$\%$"
    else:
        return s + "%"


def setup_theme(
    context="notebook",
    style="darkgrid",
    palette="deep",
    font="Helvetica",
    usetex: bool = True,
):
    try:
        extra_rc = {
            "lines.linewidth": 1,
            "lines.markeredgewidth": 1,
            "patch.edgecolor": "k",
        }
        sns.set_theme(context=context, style=style, palette=palette, rc=extra_rc)
    except (ImportError, SyntaxError):
        pass

    if usetex:
        rc("text", usetex=True)
    else:
        logger.info("Set text.usetex=%s. Font styles may be inconsistent.", usetex)
        rc("text", usetex=False)

    if font == "Helvetica":
        rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
    elif font == "Palatino":
        rc("font", **{"family": "serif", "serif": ["Palatino"]})
    elif font == "Schoolbook":
        rc("font", **{"family": "serif", "serif": ["Century Schoolbook L"]})


def asciiaxis(x, digit=1):
    if isinstance(x, int):
        x = str(x)
    elif isinstance(x, float):
        x = "{0:.{1}f}".format(x, digit)
    elif isinstance(x, np.int64):
        x = str(x)
    elif isinstance(x, np.ndarray):
        assert len(x) == 2
        x = str(x).replace("]", ")")  # upper bound not inclusive

    return x


def asciiplot(x, y, digit=1, width=50, title=None, char="="):
    """
    Print out a horizontal plot using ASCII chars.
    width is the textwidth (height) of the plot.
    """
    ax = np.array(x)
    ay = np.array(y)

    if title:
        print("[bold white]".format(title), file=sys.stderr)

    az = ay * width // ay.max()
    tx = [asciiaxis(x, digit=digit) for x in ax]
    rjust = max([len(x) for x in tx]) + 1

    for x, y, z in zip(tx, ay, az):
        x = x.rjust(rjust)
        y = y or ""
        z = "[green]{}".format(char * z)
        print("{} | {} {}".format(x, z, y), file=sys.stderr)


def print_colors(palette, outfile="Palette.png"):
    """
    print color palette (a tuple) to a PNG file for quick check
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    xmax = 20 * (len(palette) + 1)
    x1s = np.arange(0, xmax, 20)
    xintervals = [10] * len(palette)
    xx = zip(x1s, xintervals)
    ax.broken_barh(xx, (5, 10), facecolors=palette)

    ax.set_ylim(0, 20)
    ax.set_xlim(0, xmax)
    ax.set_axis_off()

    savefig(outfile)


def plot_heatmap(
    ax,
    M: np.ndarray,
    breaks: List[int],
    groups: List[Tuple[int, int, List[Tuple[int, str]], str]] = [],
    plot_breaks: bool = False,
    cmap: Optional[Union[str, Colormap]] = None,
    binsize: Optional[int] = None,
):
    """Plot heatmap illustrating the contact probabilities in Hi-C data.

    Args:
        ax (pyplot.axes): Matplotlib axis
        M (np.array): 2D numpy-array
        breaks (List[int]): Positions of chromosome starts. Can be None.
        iopts (OptionParser options): Graphical options passed in from commandline
        groups (List, optional): [(start, end, [(position, seqid)], color)]. Defaults to [].
        plot_breaks (bool): Whether to plot white breaks. Defaults to False.
        cmap (str | Colormap, optional): Colormap. Defaults to None, which uses cubehelix.
        binsize (int, optional): Resolution of the heatmap.
    """
    cmap = cmap or sns.cubehelix_palette(rot=0.5, as_cmap=True)
    ax.imshow(M, cmap=cmap, interpolation="none")
    _, xmax = ax.get_xlim()
    xlim = (0, xmax)
    if plot_breaks:
        for b in breaks[:-1]:
            ax.plot([b, b], xlim, "w-")
            ax.plot(xlim, [b, b], "w-")

    def simplify_seqid(seqid):
        seqid = seqid.replace("_", "")
        if seqid[:3].lower() == "chr":
            seqid = seqid[3:]
        return seqid.lstrip("0")

    for start, end, position_seqids, color in groups:
        # Plot a square
        ax.plot([start, start], [start, end], "-", color=color)
        ax.plot([start, end], [start, start], "-", color=color)
        ax.plot([start, end], [end, end], "-", color=color)
        ax.plot([end, end], [start, end], "-", color=color)
        for position, seqid in position_seqids:
            seqid = simplify_seqid(seqid)
            ax.text(position, end, seqid, ha="center", va="top")

    ax.set_xlim(xlim)
    ax.set_ylim((xlim[1], xlim[0]))  # Flip the y-axis so the origin is at the top
    ax.set_xticklabels(ax.get_xticks(), family="Helvetica", color="gray")
    ax.set_yticklabels(ax.get_yticks(), family="Helvetica", color="gray", rotation=90)
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    if binsize is not None:
        formatter = ticker.FuncFormatter(
            lambda x, pos: human_readable(int(x) * binsize, pos, base=True)
        )
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        title = f"Resolution = {human_size(binsize, precision=0)} per bin"
        ax.set_xlabel(title)


def discrete_rainbow(N=7, cmap=cm.Set1, usepreset=True, shuffle=False, plot=False):
    """
    Return a discrete colormap and the set of colors.

    modified from
    <http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations>

    cmap: colormap instance, eg. cm.jet.
    N: Number of colors.

    Example
    >>> x = resize(arange(100), (5,100))
    >>> djet = cmap_discretize(cm.jet, 5)
    >>> imshow(x, cmap=djet)

    See available matplotlib colormaps at:
    <http://dept.astro.lsa.umich.edu/~msshin/science/code/matplotlib_cm/>

    If N>20 the sampled colors might not be very distinctive.
    If you want to error and try anyway, set usepreset=False
    """
    import random
    from scipy import interpolate

    if usepreset:
        if 0 < N <= 5:
            cmap = cm.gist_rainbow
        elif N <= 20:
            cmap = cm.Set1
        else:
            sys.exit(discrete_rainbow.__doc__)

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = np.linspace(0, 1.0, N)
    # N+1 indices
    indices = np.linspace(0, 1.0, N + 1)
    rgbs = []
    for key in ("red", "green", "blue"):
        # Find the N colors
        D = np.array(cdict[key])
        I = interpolate.interp1d(D[:, 0], D[:, 1])
        colors = I(colors_i)
        rgbs.append(colors)
        # Place these colors at the correct indices.
        A = np.zeros((N + 1, 3), float)
        A[:, 0] = indices
        A[1:, 1] = colors
        A[:-1, 2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)

    palette = zip(*rgbs)

    if shuffle:
        random.shuffle(palette)

    if plot:
        print_colors(palette)

    # Return (colormap object, RGB tuples)
    return mpl.colors.LinearSegmentedColormap("colormap", cdict, 1024), palette


def get_intensity(octal):
    from math import sqrt

    r, g, b = octal[1:3], octal[3:5], octal[5:]
    r, g, b = int(r, 16), int(g, 16), int(b, 16)
    intensity = sqrt((r * r + g * g + b * b) / 3)
    return intensity


def adjust_spines(ax, spines, outward=False, color="lightslategray"):
    # Modified from <http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html>
    for loc, spine in ax.spines.items():
        if loc in spines:
            if outward:
                spine.set_position(("outward", 8))  # outward by 10 points
            spine.set_color(color)
        else:
            spine.set_color("none")  # don't draw spine

    if "left" in spines:
        ax.yaxis.set_ticks_position("left")
    else:
        ax.yaxis.set_ticks_position("right")

    if "bottom" in spines:
        ax.xaxis.set_ticks_position("bottom")
    else:
        ax.xaxis.set_ticks_position("top")

    # Change tick styles directly
    ax.tick_params(color=color)
    set_helvetica_axis(ax)


def set_ticklabels_helvetica(ax, xcast=int, ycast=int):
    xticklabels = [xcast(x) for x in ax.get_xticks()]
    yticklabels = [ycast(x) for x in ax.get_yticks()]
    ax.set_xticklabels(xticklabels, family="Helvetica")
    ax.set_yticklabels(yticklabels, family="Helvetica")


def draw_cmap(ax, cmap_text, vmin, vmax, cmap=None, reverse=False):
    # Draw a horizontal colormap at bottom-right corder of the canvas
    Y = np.outer(np.ones(10), np.arange(0, 1, 0.01))
    if reverse:
        Y = Y[::-1]
    xmin, xmax = 0.6, 0.9
    ymin, ymax = 0.02, 0.04
    ax.imshow(Y, extent=(xmin, xmax, ymin, ymax), cmap=cmap)
    ax.text(
        xmin - 0.01,
        (ymin + ymax) * 0.5,
        markup(cmap_text),
        ha="right",
        va="center",
        size=10,
    )
    vmiddle = (vmin + vmax) * 0.5
    xmiddle = (xmin + xmax) * 0.5
    for x, v in zip((xmin, xmiddle, xmax), (vmin, vmiddle, vmax)):
        ax.text(x, ymin - 0.005, "%.1f" % v, ha="center", va="top", size=10)


def write_messages(ax, messages, ypad=0.04):
    """
    Write text on canvas, usually on the top right corner.
    """
    tc = "gray"
    axt = ax.transAxes
    yy = 0.95
    for msg in messages:
        ax.text(0.95, yy, msg, color=tc, transform=axt, ha="right")
        yy -= ypad


def quickplot_ax(
    ax,
    data,
    xmin,
    xmax,
    xlabel,
    title=None,
    ylabel="Counts",
    counts=True,
    percentage=True,
    highlight=None,
):
    # TODO: redundant with quickplot(), need to be refactored.
    if percentage:
        total_length = sum(data.values())
        data = dict((k, v * 100.0 / total_length) for (k, v) in data.items())

    left, height = zip(*sorted(data.items()))
    pad = max(height) * 0.01
    c1, c2 = "darkslategray", "tomato"
    if counts:
        for l, h in zip(left, height):
            if xmax and l > xmax:
                break
            tag = str(int(h))
            rotation = 90
            if percentage:
                tag = append_percentage(tag) if int(tag) > 0 else ""
                rotation = 0
            color = c1
            if highlight is not None and l in highlight:
                color = c2
            ax.text(
                l,
                h + pad,
                tag,
                color=color,
                size=8,
                ha="center",
                va="bottom",
                rotation=rotation,
            )
    if xmax is None:
        xmax = max(left)

    ax.bar(left, height, align="center", fc=c1)
    if highlight:
        for h in highlight:
            ax.bar([h], [data[h]], align="center", ec=c2, fc=c2)

    ax.set_xlabel(markup(xlabel))
    if ylabel:
        ax.set_ylabel(markup(ylabel))
    if title:
        ax.set_title(markup(title))
    ax.set_xlim((xmin - 0.5, xmax + 0.5))
    if percentage:
        ax.set_ylim(0, 100)


def quickplot(
    data,
    xmin,
    xmax,
    xlabel,
    title,
    ylabel="Counts",
    figname="plot.pdf",
    counts=True,
    print_stats=True,
):
    """
    Simple plotting function - given a dictionary of data, produce a bar plot
    with the counts shown on the plot.
    """
    plt.figure(1, (6, 6))
    left, height = zip(*sorted(data.items()))
    pad = max(height) * 0.01
    if counts:
        for l, h in zip(left, height):
            if xmax and l > xmax:
                break
            plt.text(
                l,
                h + pad,
                str(h),
                color="darkslategray",
                size=8,
                ha="center",
                va="bottom",
                rotation=90,
            )
    if xmax is None:
        xmax = max(left)

    plt.bar(left, height, align="center")
    plt.xlabel(markup(xlabel))
    plt.ylabel(markup(ylabel))
    plt.title(markup(title))
    plt.xlim((xmin - 0.5, xmax + 0.5))

    # Basic statistics
    messages = []
    counts_over_xmax = sum([v for k, v in data.items() if k > xmax])
    if counts_over_xmax:
        messages += ["Counts over xmax({0}): {1}".format(xmax, counts_over_xmax)]
    kk = []
    for k, v in data.items():
        kk += [k] * v
    messages += ["Total: {0}".format(np.sum(height))]
    messages += ["Maximum: {0}".format(np.max(kk))]
    messages += ["Minimum: {0}".format(np.min(kk))]
    messages += ["Average: {0:.2f}".format(np.mean(kk))]
    messages += ["Median: {0}".format(np.median(kk))]
    ax = plt.gca()
    if print_stats:
        write_messages(ax, messages)

    set_human_axis(ax)
    set_ticklabels_helvetica(ax)
    savefig(figname)
