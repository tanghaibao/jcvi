#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import sys
import logging

from functools import partial
from glob import glob

import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib import cm, rc
from matplotlib.patches import Rectangle, Polygon, CirclePolygon, PathPatch
from matplotlib.path import Path
from matplotlib.colors import LogNorm
from matplotlib.transforms import Affine2D

from jcvi.utils.cbook import human_size
from jcvi.apps.console import dark, green


class ImageOptions (object):

    def __init__(self, opts):
        self.w, self.h = [int(x) for x in opts.figsize.split('x')]
        self.dpi = opts.dpi
        self.format = opts.format

    def __str__(self):
        return "({0}px x {1}px)".format(self.dpi * self.w, self.dpi * self.h)


def savefig(figname, dpi=150, iopts=None):
    try:
        format = figname.rsplit(".", 1)[-1].lower()
    except:
        format = "pdf"
    try:
        plt.savefig(figname, dpi=dpi, format=format)
    except:
        logging.error("savefig failed. Reset usetex to False.")
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
    return _(x)


human_readable_base = partial(human_readable, base=True)
human_formatter = ticker.FuncFormatter(human_readable)
human_base_formatter = ticker.FuncFormatter(human_readable_base)
tex_formatter = ticker.FuncFormatter(lambda x, pos: _(str(int(x))))
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


def set_image_options(instance, args=None, figsize="6x6", dpi=300,
                      format="pdf", theme="helvetica"):
    """
    Add image format options for given command line programs.
    """
    from optparse import OptionParser, OptionGroup
    assert isinstance(instance, OptionParser)

    allowed_format = ("emf", "eps", "pdf", "png", "ps", \
                      "raw", "rgba", "svg", "svgz")
    allowed_themes = ("helvetica", "palatino", "schoolbook", "mpl")

    group = OptionGroup(instance, "Image options")
    instance.add_option_group(group)

    group.add_option("--figsize", default=figsize,
            help="Figure size `width`x`height` in inches [default: %default]")
    group.add_option("--dpi", default=dpi, type="int",
            help="Physical dot density (dots per inch) [default: %default]")
    group.add_option("--format", default=format, choices=allowed_format,
            help="Generate image of format, must be one of {0}".\
            format("|".join(allowed_format)) + " [default: %default]")
    group.add_option("--theme", default=theme, choices=allowed_themes,
            help="Font theme, must be one of {0}".format("|".join(allowed_themes)) + \
                 " [default: %default]")

    if args is None:
        args = sys.argv[1:]

    opts, args = instance.parse_args(args)

    assert opts.dpi > 0
    assert "x" in opts.figsize

    setup_theme(opts.theme)

    return opts, args, ImageOptions(opts)


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
