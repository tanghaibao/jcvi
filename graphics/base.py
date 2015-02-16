#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import sys
import logging

from functools import partial

import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib import cm, rc
from matplotlib.patches import Rectangle, Polygon, CirclePolygon, PathPatch, \
            FancyArrow, FancyArrowPatch
from matplotlib.path import Path

from jcvi.utils.brewer2mpl import get_map
from jcvi.formats.base import LineFile
from jcvi.apps.console import dark, green
from jcvi.apps.base import glob, listify, datadir


class ImageOptions (object):

    def __init__(self, opts):
        self.w, self.h = [int(x) for x in opts.figsize.split('x')]
        self.dpi = opts.dpi
        self.format = opts.format
        self.cmap = cm.get_cmap(opts.cmap)
        self.opts = opts

    def __str__(self):
        return "({0}px x {1}px)".format(self.dpi * self.w, self.dpi * self.h)

    @property
    def diverge(self):
        colors = get_map(self.opts.diverge, 'diverging', 5).mpl_colors
        return colors[0], colors[-1]


class TextHandler (object):

    def __init__(self, fig):
        self.build_height_array(fig)

    @classmethod
    def get_text_width_height(cls, fig, txt="chr01", size=12, usetex=True):
        tp = mpl.textpath.TextPath((0,0), txt, size=size, usetex=usetex)
        bb = tp.get_extents()
        xmin, ymin = fig.transFigure.inverted().transform((bb.xmin, bb.ymin))
        xmax, ymax = fig.transFigure.inverted().transform((bb.xmax, bb.ymax))
        return xmax - xmin, ymax - ymin

    def build_height_array(self, fig, start=1, stop=36):
        self.heights = []
        for i in xrange(start, stop + 1):
            w, h = TextHandler.get_text_width_height(fig, size=i)
            self.heights.append((h, i))

    def select_fontsize(self, height, minsize=1, maxsize=12):
        from bisect import bisect_left

        i = bisect_left(self.heights, (height,))
        size = self.heights[i - 1][1] if i else minsize
        size = min(size, maxsize)
        return size


class AbstractLayout (LineFile):
    """
    Simple csv layout file for complex plotting settings. Typically, each line
    represents a subplot, a track or a panel.
    """
    def __init__(self, filename, delimiter=','):
        super(AbstractLayout, self).__init__(filename)

    def assign_array(self, attrib, array):
        assert len(array) == len(self)
        for x, c in zip(self, array):
            if not getattr(x, attrib):
                setattr(x, attrib, c)

    def assign_colors(self):
        colorset = get_map('Set2', 'qualitative', len(self)).mpl_colors
        self.assign_array("color", colorset)

    def assign_markers(self):
        from random import sample
        markerset = sample(mpl.lines.Line2D.filled_markers, len(self))
        self.assign_array("marker", markerset)

    def __str__(self):
        return "\n".join(str(x) for x in self)


CHARS = {
    '&':  r'\&',
    '%':  r'\%',
    '$':  r'\$',
    '#':  r'\#',
    '_':  r'\_',
    '{':  r'\{',
    '}':  r'\}',
}


def latex(s):
    return "".join([CHARS.get(char, char) for char in s])


def shorten(s, maxchar=20):
    if s <= maxchar:
        return s
    pad = (maxchar - 3) / 2
    return s[:pad] + "..." + s[-pad:]


def prettyplot():
    # Get Set2 from ColorBrewer, a set of colors deemed colorblind-safe and
    # pleasant to look at by Drs. Cynthia Brewer and Mark Harrower of Pennsylvania
    # State University. These colors look lovely together, and are less
    # saturated than those colors in Set1. For more on ColorBrewer, see:
    set2 = get_map('Set2', 'qualitative', 8).mpl_colors
    set1 = get_map('Set1', 'qualitative', 9).mpl_colors

    reds = mpl.cm.Reds
    reds.set_bad('white')
    reds.set_under('white')

    blues_r = mpl.cm.Blues_r
    blues_r.set_bad('white')
    blues_r.set_under('white')

    # Need to 'reverse' red to blue so that blue=cold=small numbers,
    # and red=hot=large numbers with '_r' suffix
    blue_red = get_map('RdBu', 'Diverging', 11,
                                  reverse=True).mpl_colormap
    green_purple = get_map('PRGn', 'diverging', 11).mpl_colormap
    red_purple = get_map('RdPu', 'Sequential', 9).mpl_colormap

    return blues_r, reds, blue_red, set1, set2, green_purple, red_purple


blues_r, reds, blue_red, set1, set2, green_purple, red_purple = prettyplot()


def normalize_axes(axes):
    axes = listify(axes)
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_axis_off()


def panel_labels(ax, labels, size=16):
    for xx, yy, panel_label in labels:
        panel_label = r"$\textbf{{{0}}}$".format(panel_label)
        ax.text(xx, yy, panel_label, size=size,
                        ha="center", va="center")


def savefig(figname, dpi=150, iopts=None):
    try:
        format = figname.rsplit(".", 1)[-1].lower()
    except:
        format = "pdf"
    try:
        plt.savefig(figname, dpi=dpi, format=format)
    except Exception as e:
        message = "savefig failed. Reset usetex to False."
        message += "\n{0}".format(str(e))
        logging.error(message)
        rc('text', **{'usetex': False})
        plt.savefig(figname, dpi=dpi)

    msg = "Figure saved to `{0}`".format(figname)
    if iopts:
        msg += " {0}".format(iopts)
    logging.debug(msg)

    plt.rcdefaults()


# human readable size (Kb, Mb, Gb)
def human_readable(x, pos, base=False):
    x = str(int(x))
    if x.endswith("000000"):
        x = x[:-6] + "M"
    elif x.endswith("000"):
        x = x[:-3] + "K"
    if base and x[-1] in "MK":
        x += "b"
    return x


human_readable_base = partial(human_readable, base=True)
human_formatter = ticker.FuncFormatter(human_readable)
human_base_formatter = ticker.FuncFormatter(human_readable_base)
mb_formatter = ticker.FuncFormatter(lambda x, pos: "{0}M".format(int(x / 1000000)))
mb_float_formatter = ticker.FuncFormatter(lambda x, pos: "{0:.1f}M".format(x / 1000000.))
kb_formatter = ticker.FuncFormatter(lambda x, pos: "{0}K".format(int(x / 1000)))


def set_human_axis(ax, formatter=human_formatter):
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)


set_human_base_axis = partial(set_human_axis, formatter=human_base_formatter)

available_fonts = [op.basename(x) for x in glob(datadir + "/*.ttf")]


def fontprop(ax, name, size=12):

    assert name in available_fonts, "Font must be one of {0}.".\
            format(available_fonts)

    import matplotlib.font_manager as fm

    fname = op.join(datadir, name)
    prop = fm.FontProperties(fname=fname, size=size)

    logging.debug("Set font to `{0}` (`{1}`).".format(name, prop.get_file()))
    for text in ax.texts:
        text.set_fontproperties(prop)

    return prop


def markup(s):
    import re
    s = re.sub("\*(.*)\*", r"\\textit{\1}", s)
    return s


def setup_theme(context='notebook', style="darkgrid", palette='deep', font='Helvetica'):
    try:
        import seaborn as sns
        extra_rc = {"lines.linewidth": 1,
                    "lines.markeredgewidth": 1,
                    "patch.edgecolor": 'k',
                    }
        sns.set(context=context, style=style, palette=palette, rc=extra_rc)
    except (ImportError, SyntaxError):
        pass

    rc('text', usetex=True)

    if font == "Helvetica":
        rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    elif font == "Palatino":
        rc('font', **{'family':'serif','serif': ['Palatino']})
    elif font == "Schoolbook":
        rc('font', **{'family':'serif','serif': ['Century Schoolbook L']})


setup_theme()


def asciiaxis(x, digit=1):
    if isinstance(x, int):
        x = str(x)
    elif isinstance(x, float):
        x = "{0:.{1}f}".format(x, digit)
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
        print >> sys.stderr, dark(title)

    az = ay * width / ay.max()
    tx = [asciiaxis(x, digit=digit) for x in ax]
    rjust = max([len(x) for x in tx]) + 1

    for x, y, z in zip(tx, ay, az):
        x = x.rjust(rjust)
        y = y or ""
        z = green(char * z)
        print >> sys.stderr, "{0} |{1} {2}".format(x, z, y)


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


def discrete_rainbow(N=7, cmap=cm.Set1, usepreset=True, shuffle=False, \
                     plot=False):
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
    colors_i = np.linspace(0,1.,N)
    # N+1 indices
    indices = np.linspace(0,1.,N+1)
    rgbs = []
    for key in ('red','green','blue'):
       # Find the N colors
       D = np.array(cdict[key])
       I = interpolate.interp1d(D[:,0], D[:,1])
       colors = I(colors_i)
       rgbs.append(colors)
       # Place these colors at the correct indices.
       A = np.zeros((N+1,3), float)
       A[:,0] = indices
       A[1:,1] = colors
       A[:-1,2] = colors
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
    return mpl.colors.LinearSegmentedColormap('colormap',cdict,1024), palette


def get_intensity(octal):
    from math import sqrt

    r, g, b = octal[1:3], octal[3:5], octal[5:]
    r, g, b = int(r, 16), int(g, 16), int(b, 16)
    intensity = sqrt((r * r + g * g + b * b) / 3)
    return intensity


def adjust_spines(ax, spines):
    # Modified from <http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html>
    for loc, spine in ax.spines.items():
        if loc in spines:
            pass
            #spine.set_position(('outward', 10)) # outward by 10 points
            #spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine

    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks_position('right')

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks_position('top')


def set_ticklabels_helvetica(ax, xcast=int, ycast=int):
    xticklabels = [xcast(x) for x in ax.get_xticks()]
    yticklabels = [ycast(x) for x in ax.get_yticks()]
    ax.set_xticklabels(xticklabels, family='Helvetica')
    ax.set_yticklabels(yticklabels, family='Helvetica')


def draw_cmap(ax, cmap_text, vmin, vmax, cmap=None, reverse=False):
    # Draw a horizontal colormap at bottom-right corder of the canvas
    Y = np.outer(np.ones(10), np.arange(0, 1, 0.01))
    if reverse:
        Y = Y[::-1]
    xmin, xmax = .6, .9
    ymin, ymax = .02, .04
    ax.imshow(Y, extent=(xmin, xmax, ymin, ymax), cmap=cmap)
    ax.text(xmin - .01, (ymin + ymax) * .5, cmap_text,
            ha="right", va="center", size=10)
    vmiddle = (vmin + vmax) * .5
    xmiddle = (xmin + xmax) * .5
    for x, v in zip((xmin, xmiddle, xmax), (vmin, vmiddle, vmax)):
        ax.text(x, ymin - .005, "%.1f" % v, ha="center", va="top", size=10)
