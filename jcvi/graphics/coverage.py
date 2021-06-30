#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog chrC01 chr.sizes data

Read coverage histogram, similar to wiggle plot. Data contains all the track
data in the form of tab-delimited (x, y) lists.
"""

import os.path as op
import sys
import logging

import numpy as np

from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import (
    plt,
    savefig,
    Rectangle,
    mb_formatter,
    mb_float_formatter,
    adjust_spines,
    get_map,
)
from jcvi.apps.base import OptionParser, glob


class XYtrack(object):
    def __init__(self, ax, datafile, color=None, ymax=40):
        self.ax = ax
        self.xy = []
        fp = open(datafile)
        for row in fp:
            atoms = row.split()
            self.xy.append([int(atoms[0]), float(atoms[1])])
        fp.close()

        self.x, self.y = zip(*self.xy)
        logging.debug(
            "File `{0}` imported (records={1}).".format(datafile, len(self.x))
        )
        self.color = color or "k"
        self.ymax = ymax

    @property
    def mapping(self):
        return dict(zip(self.x, self.y))

    def interpolate(self, maxsize, unit=10000):
        maxsize = int(maxsize)
        for pos in range(unit, maxsize + unit, unit):
            if pos in self.x:
                continue
            self.xy.append([pos, 0])
        self.xy.sort()
        self.x, self.y = zip(*self.xy)
        logging.debug("After interpolate: {0}".format(len(self.x)))

    def cap(self, ymax):
        self.xy = [[a, 0] if b > ymax else [a, b] for a, b in self.xy]
        self.x, self.y = zip(*self.xy)

    def draw(self):
        ax = self.ax
        color = self.color
        ax.plot(self.x, self.y, lw=0)
        ax.fill_between(self.x, self.y, color=color, lw=0)
        ax.set_ylim(0, self.ymax)
        ax.set_axis_off()

    def import_hlfile(self, hlfile, chr, unit=10000, diverge=("r", "g")):
        rr, gg = diverge
        fp = open(hlfile)
        imported = 0
        mapping = self.mapping
        for row in fp:
            if row.strip() == "":
                continue
            seqid, start, end, tag = row.split()
            if seqid != chr:
                continue
            start = int(start)
            end = int(end)
            if tag == "double":
                self.highlight(mapping, start, end, color=rr, unit=unit)
            else:
                self.highlight(mapping, start, end, color=gg, unit=unit)
            imported += 1
        logging.debug("Imported {0} regions from file `{1}`.".format(imported, hlfile))

    def highlight(self, mapping, start, end, color="r", unit=10000, zorder=10):
        ax = self.ax
        x = range(start, end + unit, unit)
        y = [mapping[z] for z in x]
        # Mask the highlight region so that they don't appear in background
        for a in self.xy:
            if start <= a[0] <= end:
                a[1] = 0
        self.x, self.y = zip(*self.xy)
        ax.plot(x, y, lw=0)
        ax.fill_between(x, y, color=color, lw=0, zorder=zorder)

    def vlines(self, xs, color="m"):
        for x in xs:
            self.ax.plot((x, x), (0, self.ymax), "-", color=color, lw=2)


class Coverage(object):
    def __init__(
        self,
        fig,
        root,
        canvas,
        chr,
        xlim,
        datadir,
        order=None,
        hlsuffix=None,
        palette=None,
        cap=50,
        gauge="bottom",
        plot_label=True,
        plot_chr_label=True,
        gauge_step=5000000,
        vlines=None,
        labels_dict={},
        diverge=("r", "g"),
    ):
        x, y, w, h = canvas
        p = 0.01
        root.add_patch(
            Rectangle(
                (x - p, y - p),
                w + 2 * p,
                h + 2 * p,
                lw=1,
                fill=False,
                ec="darkslategray",
                zorder=10,
            )
        )
        datafiles = glob(op.join(datadir, chr + "*"))

        if order:
            datafiles = [z for z in datafiles if z.split(".")[1] in order]
            datafiles.sort(key=lambda x: order.index(x.split(".")[1]))

        ntracks = len(datafiles)
        yinterval = h / ntracks
        yy = y + h

        if palette is None:
            # Get the palette
            set2 = get_map("Set2", "qualitative", ntracks).mpl_colors
        else:
            set2 = [palette] * ntracks

        if gauge == "top":
            gauge_ax = fig.add_axes([x, yy + p, w, 0.0001])
            adjust_spines(gauge_ax, ["top"])
            tpos = yy + 0.07
        elif gauge == "bottom":
            gauge_ax = fig.add_axes([x, y - p, w, 0.0001])
            adjust_spines(gauge_ax, ["bottom"])
            tpos = y - 0.07

        start, end = xlim
        if gauge:
            fs = gauge_step < 1000000
            setup_gauge_ax(gauge_ax, start, end, gauge_step, float_formatter=fs)

        if plot_chr_label:
            root.text(
                x + w / 2,
                tpos,
                chr,
                ha="center",
                va="center",
                color="darkslategray",
                size=16,
            )

        yys = []
        for label, datafile, c in zip(order, datafiles, set2):
            yy -= yinterval
            yys.append(yy)
            ax = fig.add_axes([x, yy, w, yinterval * 0.9])
            xy = XYtrack(ax, datafile, color=c)
            xy.interpolate(end)
            xy.cap(ymax=cap)
            if vlines:
                xy.vlines(vlines)
            if hlsuffix:
                hlfile = op.join(datadir, ".".join((label, hlsuffix)))
                xy.import_hlfile(hlfile, chr, diverge=diverge)
            if plot_label:
                label = labels_dict.get(label, label.capitalize())
                label = r"\textit{{{0}}}".format(label)
                root.text(x - 0.015, yy + yinterval / 2, label, ha="right", va="center")
            xy.draw()
            ax.set_xlim(*xlim)

        self.yys = yys


def setup_gauge_ax(gauge_ax, start, end, gauge_step, float_formatter=False):
    gauge_ax.set_xlim(start, end)
    formatter = mb_float_formatter if float_formatter else mb_formatter
    gauge_ax.xaxis.set_major_formatter(formatter)
    gauge_ax.xaxis.set_ticks(np.arange(start + gauge_step, end, gauge_step))
    gauge_ax.yaxis.set_ticks([])


def main():
    p = OptionParser(__doc__)
    p.add_option("--order", help="The order to plot the tracks, comma-separated")
    opts, args, iopts = p.set_image_options()

    if len(args) != 3:
        sys.exit(not p.print_help())

    chr, sizes, datadir = args
    order = opts.order
    hlsuffix = opts.hlsuffix
    if order:
        order = order.split(",")
    sizes = Sizes(sizes)
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    canvas = (0.12, 0.35, 0.8, 0.35)
    chr_size = sizes.get_size(chr)
    Coverage(
        fig, root, canvas, chr, (0, chr_size), datadir, order=order, hlsuffix=hlsuffix
    )

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = chr + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main()
