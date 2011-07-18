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


def asciiplot(x, y, width=50, title=None, char="="):
    """
    Print out a horizontal plot using ASCII chars.
    width is the textwidth (height) of the plot.
    """
    ax = np.array(x)
    ay = np.array(y)

    if title:
        print >> sys.stderr, ColoredText(title, "dark")

    az = ay * width / ay.max()
    for x, y, z in zip(ax, ay, az):
        x = str(x) if isinstance(x, int) else "{0:.1f}".format(x)
        x = x.rjust(10)
        y = y or ""
        z = ColoredText(char * z, "green")
        print >> sys.stderr, "{0} |{1} {2}".format(x, z, y)
