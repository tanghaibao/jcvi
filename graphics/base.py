#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys

from functools import partial

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle, Polygon, CirclePolygon
from matplotlib.colors import LogNorm
from matplotlib import cm

from jcvi.utils.cbook import human_size
from jcvi.apps.console import dark, green

# i always like the latex font
_ = lambda x: r"$\mathsf{%s}$" % str(x).replace("_", "\_").replace(" ", r"\ ")


class ImageOptions (object):

    def __init__(self, opts):
        self.w, self.h = [int(x) for x in opts.figsize.split('x')]
        self.dpi = opts.dpi
        self.format = opts.format

    def __str__(self):
        return "({0}px x {1}px)".format(self.dpi * self.w, self.dpi * self.h)


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


def set_tex_axis(ax, formatter=tex_formatter):
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)


set_human_axis = partial(set_tex_axis, formatter=human_formatter)
set_human_base_axis = partial(set_tex_axis, formatter=human_base_formatter)


def set_image_options(instance, args=None, figsize="6x6", dpi=300, format="pdf"):
    """
    Add image format options for given command line programs.
    """
    from optparse import OptionParser, OptionGroup
    assert isinstance(instance, OptionParser)

    allowed_format = ("emf", "eps", "pdf", "png", "ps", \
                      "raw", "rgba", "svg", "svgz")

    group = OptionGroup(instance, "Image options")
    instance.add_option_group(group)

    group.add_option("--figsize", default=figsize,
            help="Figure size `width`x`height` in inches [default: %default]")
    group.add_option("--dpi", default=dpi, type="int",
            help="Physical dot density (dots per inch) [default: %default]")
    group.add_option("--format", default=format, choices=allowed_format,
            help="Generate image of format, must be one of {0}".\
            format("|".join(allowed_format)) + " [default: %default]")

    if args is None:
        args = sys.argv[1:]

    opts, args = instance.parse_args(args)

    assert opts.dpi > 0
    assert "x" in opts.figsize

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
