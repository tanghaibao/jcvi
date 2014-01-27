#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog napuscov chrC01 chr.sizes data

Read coverage histogram, similar to wiggle plot. Data contains all the track
data in the form of tab-delimited (x, y) lists.
"""

import os.path as op
import sys
import logging

import numpy as np

from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import plt, savefig, Rectangle, mb_formatter, \
        adjust_spines
from jcvi.apps.base import OptionParser, debug, glob
debug()


class XYtrack (object):

    def __init__(self, datafile, color=None):
        self.x, self.y = np.loadtxt(datafile, unpack=True)
        logging.debug("File `{0}` imported (records={1})."\
                        .format(datafile, len(self.x)))
        self.color = color

    def draw(self, ax, label):
        x, y = self.x, self.y
        color = self.color
        ax.plot(x, y, color=color)
        ax.fill_between(x, y, color=color)
        ax.set_ylim(0, 40)
        ax.set_axis_off()


class Coverage (object):

    def __init__(self, fig, root, canvas, chr, xlim, datadir, order=None,
                 gauge="bottom"):
        x, y, w, h = canvas
        p = .01
        root.add_patch(Rectangle((x - p, y - p), w + 2 * p, h + 2 * p, lw=1,
                        fill=False, ec="darkslategray", zorder=10))
        datafiles = glob(op.join(datadir, chr + "*"))
        ntracks = len(datafiles)
        yinterval = h / ntracks
        yy = y + h

        # Get the pallette
        try:
            import brewer2mpl
            set2 = brewer2mpl.get_map('Set2', 'qualitative', ntracks).mpl_colors
        except ImportError:
            from jcvi.graphics.base import set2

        if order:
            datafiles.sort(key=lambda x: order.index(x.split(".")[1]))

        if gauge == "top":
            gauge_ax = fig.add_axes([x, yy + p, w, .0001])
            adjust_spines(gauge_ax, ["top"])
            tpos = yy + .07
        elif gauge == "bottom":
            gauge_ax = fig.add_axes([x, y - p, w, .0001])
            adjust_spines(gauge_ax, ["bottom"])
            tpos = y - .07

        gauge_ax.set_xlim(*xlim)
        gauge_ax.xaxis.set_major_formatter(mb_formatter)
        gauge_ax.yaxis.set_ticks([])

        root.text(x + w / 2, tpos, chr, ha="center", va="center",
                  color="darkslategray")

        for datafile, c in zip(datafiles, set2):
            yy -= yinterval
            ax = fig.add_axes([x, yy, w, yinterval * .9])
            xy = XYtrack(datafile, color=c)
            label = datafile.split(".")[1]
            xy.draw(ax, label)
            ax.set_xlim(*xlim)
            root.text(x - .05, yy + yinterval / 2, label,
                        ha="center", va="center", color=c)


def main():
    p = OptionParser(__doc__)
    p.add_option("--order",
                help="The order to plot the tracks, comma-separated")
    opts, args, iopts = p.set_image_options()

    if len(args) != 3:
        sys.exit(not p.print_help())

    chr, sizes, datadir = args
    order = opt.order
    if order:
        order = order.split(",")
    sizes = Sizes(sizes)
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    canvas = (.12, .35, .8, .35)
    chr_size = sizes.get_size(chr)
    c = Coverage(fig, root, canvas, chr, (0, chr_size), datadir,
                 order=order)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = chr + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
