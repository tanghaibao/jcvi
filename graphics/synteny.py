#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog mcscan.txt all.bed layout.csv

Illustrate MCscan multiple collinearity alignments. Use layout.csv to indicate
the positions of tracks. For example:

#x, y, rotation, label_position, color
0.5, 0.6, 0, left, g
0.25, 0.7, 45, top, m

With the row ordering corresponding to the column ordering in the MCscan output.
"""


import sys
import logging
import numpy as np

from optparse import OptionParser

from jcvi.algorithms.synteny import BlockFile
from jcvi.formats.bed import Bed
from jcvi.formats.base import LineFile
from jcvi.graphics.base import plt, _, set_image_options, Affine2D
from jcvi.utils.cbook import human_size
from jcvi.apps.base import debug
debug()


class LayoutLine (object):

    def __init__(self, row, delimiter=','):
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.x = float(args[0])
        self.y = float(args[1])
        self.rotation = int(args[2])
        self.label_position = args[3]
        self.color = args[4]


class Layout (LineFile):

    def __init__(self, filename, delimiter=','):
        super(Layout, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            if row[0] == '#':
                continue
            self.append(LayoutLine(row, delimiter=delimiter))


class Region (object):

    def __init__(self, ax, ext, layout, pad=.05, aspect_ratio=1):
        x, y = layout.x, layout.y
        lr = layout.rotation
        tr = Affine2D().translate(-x, -y).rotate_deg(lr).\
                        translate(x, y) + ax.transAxes

        ax.plot((x - .15, x + .15), (y, y), color="gray", transform=tr, lw=2)

        start, end, si, ei, chr, orientation = ext
        startbp, endbp = human_size(start.start), human_size(end.end)
        if orientation == '-':
            startbp, endbp = endbp, startbp
        label = "{0}:{1}-{2}".format(chr, startbp, endbp)

        va, ha = "center", "center"
        lp = layout.label_position
        if lp in ("left", "right"):
            ha = lp
        else:
            va = lp

        l = np.array((x, y + pad))
        trans_angle = ax.transAxes.transform_angles(np.array((lr, )),
                                                    l.reshape((1, 2)))[0]

        t = ax.text(l[0], l[1], _(label), color=layout.color,
                    ha=ha, va=va, size=10, rotation=trans_angle)


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = set_image_options(p, figsize="8x6")

    if len(args) != 3:
        sys.exit(not p.print_help())

    datafile, bedfile, layoutfile = args
    bed = Bed(bedfile)
    order = bed.order
    bf = BlockFile(datafile)
    lo = Layout(layoutfile)

    pf = datafile.rsplit(".", 1)[0]
    aspect_ratio = iopts.w * 1. / iopts.h
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    for i in xrange(bf.ncols):
        ext = bf.get_extent(i, order)
        Region(root, ext, lo[i], aspect_ratio=aspect_ratio)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


if __name__ == '__main__':
    main()
