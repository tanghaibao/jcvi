#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from matplotlib import cm

from jcvi.utils.cbook import human_size

# i always like the latex font
_ = lambda x: r"$\mathsf{%s}$" % x.replace("_", " ").replace(" ", r"\ ")

# human readable size (Kb, Mb, Gb)
human_size_formatter = ticker.FuncFormatter(lambda x, pos: human_size(x, False,
    precision=0))


