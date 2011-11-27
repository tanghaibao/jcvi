#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog csvfile

Draw heatmap based on the data in the csv file. In a microarray setting, the
rows represent genes, and columns represent conditions. Some conditions can be
grouped which the script expect to see on the first row when --groups is on.
Something like::

,WT+BL,,,,irx8+BL,,,,OE+BL,,,,WT,,,,irx8,,,,OE,,,
, Day 0,Day 3,Day 6,Day 9, Day 0,Day 3,Day 6,Day 9, Day 0,Day 3,Day 6,Day 9, ...
GAUT12,0.801069878,15.34822591,5.897076869,26.17286587,0,0,0,0,296.1121751, ...
MYB46,0.812252396,31.12495832,11.39240156,44.63179732,4.469148552,57.28160454, ...
"""


import sys
import logging

import numpy as np
from itertools import groupby
from optparse import OptionParser

from jcvi.graphics.base import plt, cm, _, set_image_options, LogNorm
from jcvi.apps.base import debug
debug()


def parse_csv(csvfile, vmin=0, groups=False):
    import csv

    reader = csv.reader(open(csvfile))
    if groups:
        groups = reader.next()[1:]
        # Fill in empty cells in groups
        filled_groups = []
        lastg = ""
        for g in groups:
            g = g.strip() or lastg
            filled_groups.append(g)
            lastg = g
        groups = filled_groups

    rows = []
    cols = reader.next()[1:]
    data = []
    for row in reader:
        name = row[0]
        d = [max(vmin, float(x)) for x in row[1:]]
        rows.append(name)
        data.append(d)

    data = np.array(data)

    return groups, rows, cols, data


def main():
    p = OptionParser(__doc__)
    p.add_option("--groups", default=False, action="store_true",
                 help="The first row contains group info [default: %default]")
    p.add_option("--cmap", default="jet",
                 help="Use this color map [default: %default]")
    opts, args, iopts = set_image_options(p, figsize="8x8")

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]

    groups, rows, cols, data = parse_csv(datafile, vmin=1, groups=opts.groups)
    cols = [x.replace("ay ", "") for x in cols]

    plt.rcParams["axes.linewidth"] = 0

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    ax = fig.add_axes([.15, .15, .7, .7])

    default_cm = cm.get_cmap(opts.cmap)
    im = ax.matshow(data, cmap=default_cm, norm=LogNorm(vmin=1, vmax=10000))
    nrows, ncols = len(rows), len(cols)

    plt.xticks(range(ncols), cols, rotation=45, size=10, ha="center")
    plt.yticks(range(nrows), rows, size=10)

    for x in ax.get_xticklines() + ax.get_yticklines():
        x.set_visible(False)

    ax.set_xlim(-.5, ncols - .5)

    t = [1, 10, 100, 1000, 10000]
    axcolor = fig.add_axes([.9, .3, .02, .4])
    fig.colorbar(im, cax=axcolor, ticks=t, format=_("%d"))

    if groups:
        groups = [(key, len(list(nn))) for key, nn in groupby(groups)]
        xstart = .15
        yy = min(.83, .5 + .5 * nrows / ncols * .7 + .06)
        xinterval = .7 / ncols
        e = .005
        sep = -.5

        for k, kl in groups:
            # Separator in the array area
            sep += kl
            ax.plot([sep, sep], [-.5, nrows - .5], "w-", lw=2)
            # Group labels on the top
            kl *= xinterval
            root.plot([xstart + e, xstart + kl - e], [yy, yy], "-", color="gray", lw=2)
            root.text(xstart + .5 * kl, yy + e, k, ha="center", color="gray")
            xstart += kl

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + opts.cmap + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


if __name__ == '__main__':
    main()
