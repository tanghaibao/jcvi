#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from functools import partial

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from matplotlib import cm

from jcvi.utils.cbook import human_size

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
