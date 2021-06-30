#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog csvfile

Draw heatmap based on the data in the csv file. In a microarray setting, the
rows represent genes, and columns represent conditions. Some conditions can be
grouped which the script expect to see on the first row when --groups is on::

,WT+BL,,,,irx8+BL,,,,OE+BL,,,,WT,,,,irx8,,,,OE,,,
, Day 0,Day 3,Day 6,Day 9, Day 0,Day 3,Day 6,Day 9, Day 0,Day 3,Day 6,Day 9, ...
GAUT12,0.801069878,15.34822591,5.897076869,26.17286587,0,0,0,0,296.1121751, ...
MYB46,0.812252396,31.12495832,11.39240156,44.63179732,4.469148552,57.28160454, ...

Option --rowgroups requires an additional file that group the genes::

I	MYB46,GUX1
II	I14H/IRX14-L,IRX10
III	I9H/IRX9-L,IRX14
IV	IRX7,GUX2
"""


import sys

import numpy as np
from itertools import groupby

from jcvi.graphics.base import mpl, plt, savefig
from jcvi.apps.base import OptionParser


def parse_csv(csvfile, vmin=0, groups=False):
    import csv

    reader = csv.reader(open(csvfile))
    if groups:
        groups = next(reader)[1:]
        # Fill in empty cells in groups
        filled_groups = []
        lastg = ""
        for g in groups:
            g = g.strip() or lastg
            filled_groups.append(g)
            lastg = g
        groups = filled_groups

    rows = []
    cols = next(reader)[1:]
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
    p.add_option(
        "--groups",
        default=False,
        action="store_true",
        help="The first row contains group info",
    )
    p.add_option("--rowgroups", help="Row groupings")
    p.add_option(
        "--horizontalbar",
        default=False,
        action="store_true",
        help="Horizontal color bar [default: vertical]",
    )
    opts, args, iopts = p.set_image_options(figsize="8x8")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (datafile,) = args
    pf = datafile.rsplit(".", 1)[0]
    rowgroups = opts.rowgroups

    groups, rows, cols, data = parse_csv(datafile, vmin=1, groups=opts.groups)
    cols = [x.replace("ay ", "") for x in cols]

    if rowgroups:
        fp = open(rowgroups)
        rgroups = []
        for row in fp:
            a, b = row.split()
            irows = [rows.index(x) for x in b.split(",")]
            rgroups.append((a, min(irows), max(irows)))

    plt.rcParams["axes.linewidth"] = 0

    xstart = 0.18
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    ax = fig.add_axes([xstart, 0.15, 0.7, 0.7])

    im = ax.matshow(data, cmap=iopts.cmap, norm=mpl.colors.LogNorm(vmin=1, vmax=10000))
    nrows, ncols = len(rows), len(cols)

    xinterval = 0.7 / ncols
    yinterval = 0.7 / max(nrows, ncols)

    plt.xticks(range(ncols), cols, rotation=45, size=10, ha="center")
    plt.yticks(range(nrows), rows, size=10)

    for x in ax.get_xticklines() + ax.get_yticklines():
        x.set_visible(False)

    ax.set_xlim(-0.5, ncols - 0.5)

    t = [1, 10, 100, 1000, 10000]
    pad = 0.06
    if opts.horizontalbar:
        ypos = 0.5 * (1 - nrows * yinterval) - pad
        axcolor = fig.add_axes([0.3, ypos, 0.4, 0.02])
        orientation = "horizontal"
    else:
        axcolor = fig.add_axes([0.9, 0.3, 0.02, 0.4])
        orientation = "vertical"
    fig.colorbar(im, cax=axcolor, ticks=t, orientation=orientation)

    if groups:
        groups = [(key, len(list(nn))) for key, nn in groupby(groups)]
        yy = 0.5 + 0.5 * nrows / ncols * 0.7 + 0.06
        e = 0.005
        sep = -0.5

        for k, kl in groups:
            # Separator in the array area
            sep += kl
            ax.plot([sep, sep], [-0.5, nrows - 0.5], "w-", lw=2)
            # Group labels on the top
            kl *= xinterval
            root.plot([xstart + e, xstart + kl - e], [yy, yy], "-", color="gray", lw=2)
            root.text(xstart + 0.5 * kl, yy + e, k, ha="center", color="gray")
            xstart += kl

    if rowgroups:
        from jcvi.graphics.glyph import TextCircle

        xpos = 0.04
        tip = 0.015
        assert rgroups
        ystart = 1 - 0.5 * (1 - nrows * yinterval)
        for gname, start, end in rgroups:
            start = ystart - start * yinterval
            end = ystart - (end + 1) * yinterval
            start -= tip / 3
            end += tip / 3

            # Bracket the groups
            root.plot((xpos, xpos + tip), (start, start), "k-", lw=2)
            root.plot((xpos, xpos), (start, end), "k-", lw=2)
            root.plot((xpos, xpos + tip), (end, end), "k-", lw=2)
            TextCircle(root, xpos, 0.5 * (start + end), gname)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + opts.cmap + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main()
