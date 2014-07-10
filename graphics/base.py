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
            FancyArrowPatch
from matplotlib.path import Path
from matplotlib.colors import LogNorm
from matplotlib.transforms import Affine2D

from jcvi.utils.brewer2mpl import get_map
from jcvi.apps.console import dark, green
from jcvi.apps.base import glob, listify


class ImageOptions (object):

    def __init__(self, opts):
        self.w, self.h = [int(x) for x in opts.figsize.split('x')]
        self.dpi = opts.dpi
        self.format = opts.format

    def __str__(self):
        return "({0}px x {1}px)".format(self.dpi * self.w, self.dpi * self.h)


def diverge_colors(scheme, color_class=5):
    colors = get_map(scheme, 'diverging', color_class).mpl_colors
    return colors[0], colors[-1]


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


def prettyplot():
    # Get Set2 from ColorBrewer, a set of colors deemed colorblind-safe and
    # pleasant to look at by Drs. Cynthia Brewer and Mark Harrower of Pennsylvania
    # State University. These colors look lovely together, and are less
    # saturated than those colors in Set1. For more on ColorBrewer, see:
    # - Flash-based interactive map:
    #     http://colorbrewer2.org/
    # - A quick visual reference to every ColorBrewer scale:
    #     http://bl.ocks.org/mbostock/5577023
    set2 = get_map('Set2', 'qualitative', 8).mpl_colors

    # Another ColorBrewer scale. This one has nice "traditional" colors like
    # reds and blues
    set1 = get_map('Set1', 'qualitative', 9).mpl_colors

    # Set some commonly used colors
    almost_black = '#262626'
    light_grey = np.array([float(248) / float(255)] * 3)

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

    # Default "patches" like scatterplots
    mpl.rcParams['patch.linewidth'] = 0.75     # edge width in points

    # Default empty circle with a colored outline
    mpl.rcParams['patch.facecolor'] = 'none'
    mpl.rcParams['patch.edgecolor'] = set2[0]

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    mpl.rcParams['axes.edgecolor'] = almost_black
    mpl.rcParams['axes.labelcolor'] = almost_black
    mpl.rcParams['axes.linewidth'] = 0.5

    # Make the default grid be white so it "removes" lines rather than adds
    mpl.rcParams['grid.color'] = 'white'

    # change the tick colors also to the almost black
    mpl.rcParams['ytick.color'] = almost_black
    mpl.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    mpl.rcParams['text.color'] = almost_black
    return almost_black, blues_r, reds, blue_red, set1, set2, \
            light_grey, green_purple, red_purple


# Code borrowed from https://github.com/olgabot/prettyplotlib (thanks)
almost_black, blues_r, reds, blue_red, set1, set2, \
    light_grey, green_purple, red_purple = prettyplot()


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
tex_formatter = ticker.FuncFormatter(lambda x, pos: _(str(int(x))))
mb_formatter = ticker.FuncFormatter(lambda x, pos: "{0}M".format(int(x / 1000000)))
mb_float_formatter = ticker.FuncFormatter(lambda x, pos: "{0:.1f}M".format(x / 1000000.))
kb_formatter = ticker.FuncFormatter(lambda x, pos: "{0}K".format(int(x / 1000)))
tex_1digit_formatter = ticker.FuncFormatter(lambda x, pos: _("{0:.1f}".format(x)))
tex_2digit_formatter = ticker.FuncFormatter(lambda x, pos: _("{0:.2f}".format(x)))


def set_tex_axis(ax, formatter=tex_formatter):
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)


set_human_axis = partial(set_tex_axis, formatter=human_formatter)
set_human_base_axis = partial(set_tex_axis, formatter=human_base_formatter)

font_dir = op.join(op.dirname(__file__), "fonts")
available_fonts = [op.basename(x) for x in glob(font_dir + "/*")]


def fontprop(ax, name, size=12):

    assert name in available_fonts, "Font must be one of {0}.".\
            format(available_fonts)

    import matplotlib.font_manager as fm

    fname = op.join(font_dir, name)
    prop = fm.FontProperties(fname=fname, size=size)

    logging.debug("Set font to `{0}` (`{1}`).".format(name, prop.get_file()))
    for text in ax.texts:
        text.set_fontproperties(prop)

    return prop


def markup(s):
    import re
    s = re.sub("\*(.*)\*", r"\\textit{\1}", s)
    return s


def setup_theme(theme="helvetica"):

    plt.rcdefaults()
    # i always like the latex font
    _ = lambda m: "\n".join(r"$\mathsf{%s}$" % str(x).\
              replace("_", "\_").replace(" ", r"\ ") for x in m.split("\n"))

    if theme == "mpl":
        return _

    rc('text', **{'usetex': True})

    if theme == "helvetica":
        rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    elif theme == "palatino":
        rc('font', **{'family':'serif','serif': ['Palatino']})
    elif theme == "schoolbook":
        rc('font', **{'family':'serif','serif': ['Century Schoolbook L']})

    return str


_ = setup_theme()


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


def cmap_map(function, cmap):
    """
    Recipe taken from:
    <http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations>

    Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous
    points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}

    # First get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))

    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step: np.array(cmap(step)[0:3])
    old_LUT = np.array(map(reduced_cmap, step_list))
    new_LUT = np.array(map(function, old_LUT))

    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red', 'green', 'blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j, i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]

        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector

    return mpl.colors.LinearSegmentedColormap('colormap', cdict, 1024)


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
