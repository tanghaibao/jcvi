#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys

from functools import partial

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle, Polygon, CirclePolygon
from matplotlib import cm

from jcvi.utils.cbook import human_size
from jcvi.apps.console import ColoredText

# i always like the latex font
_ = lambda x: r"$\mathsf{%s}$" % str(x).replace("_", " ").replace(" ", r"\ ")

# human readable size (Kb, Mb, Gb)
human_size_formatter = ticker.FuncFormatter(lambda x, pos: \
        _(human_size(x, precision=0)))
tex_formatter = ticker.FuncFormatter(lambda x, pos: _(str(int(x))))


def set_tex_axis(ax, formatter=tex_formatter):
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

set_human_axis = partial(set_tex_axis, formatter=human_size_formatter)


def set_format(instance, default="pdf"):
    """
    Add image format options for given command line programs.
    """
    from optparse import OptionParser
    assert isinstance(instance, OptionParser)

    allowed_format = ("emf", "eps", "pdf", "png", "ps", \
                      "raw", "rgba", "svg", "svgz")

    instance.add_option("--format", default=default,
            help="Generate image of format, must be one of {0}".\
            format("|".join(allowed_format)) + " [default: %default]")


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
        print >> sys.stderr, ColoredText(title, "dark")

    az = ay * width / ay.max()
    tx = [asciiaxis(x, digit=digit) for x in ax]
    rjust = max([len(x) for x in tx]) + 1

    for x, y, z in zip(tx, ay, az):
        x = x.rjust(rjust)
        y = y or ""
        z = ColoredText(char * z, "green")
        print >> sys.stderr, "{0} |{1} {2}".format(x, z, y)
