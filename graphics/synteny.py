#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog mcscan.txt all.bed layout.csv

Illustrate MCscan multiple collinearity alignments. Use layout.csv to indicate
the positions of tracks. For example:

#x, y, rotation, ha, va, color
0.5, 0.6, 0, left, center, g
0.25, 0.7, 45, top, center, m

With the row ordering corresponding to the column ordering in the MCscan output.
"""


import sys
import logging
import numpy as np

from optparse import OptionParser

from jcvi.algorithms.synteny import BlockFile
from jcvi.formats.bed import Bed
from jcvi.formats.base import LineFile, DictFile
from jcvi.graphics.base import plt, _, set_image_options, Affine2D, \
        Path, PathPatch
from jcvi.graphics.glyph import Glyph
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
        self.ha = args[3]
        self.va = args[4]
        self.color = args[5]


class Layout (LineFile):

    def __init__(self, filename, delimiter=','):
        super(Layout, self).__init__(filename)
        fp = open(filename)
        self.edges = []
        for row in fp:
            if row[0] == '#':
                continue
            if row[0] == 'e':
                args = row.rstrip().split(delimiter)
                args = [x.strip() for x in args]
                a, b = args[1:3]
                a, b = int(a), int(b)
                assert args[0] == 'e'
                self.edges.append((a, b))
            else:
                self.append(LayoutLine(row, delimiter=delimiter))


class Shade (object):

    def __init__(self, ax, a, b, ymid, highlight=False):
        a1, a2 = a
        b1, b2 = b
        ax1, ay1 = a1
        ax2, ay2 = a2
        bx1, by1 = b1
        bx2, by2 = b2
        M, C4, L, CP = Path.MOVETO, Path.CURVE4, Path.LINETO, Path.CLOSEPOLY
        pathdata = \
        [
            (M, a1),
            (C4, (ax1, ymid)),
            (C4, (bx1, ymid)),
            (C4, b1),
            (L, b2),
            (C4, (bx2, ymid)),
            (C4, (ax2, ymid)),
            (C4, a2),
            (CP, a1)
        ]
        codes, verts = zip(*pathdata)
        path = Path(verts, codes)
        if highlight:
            ec = fc = 'orange'
            alpha = .8
        else:
            ec = fc = "k"
            alpha = .2

        ax.add_patch(PathPatch(path, fc=fc, ec=ec, alpha=alpha))


class Region (object):

    def __init__(self, ax, ext, layout, bed, scale, switch=None, pad=.04):
        x, y = layout.x, layout.y
        self.y = y
        lr = layout.rotation
        tr = Affine2D().rotate_deg_around(x, y, lr) + ax.transAxes
        inv = ax.transAxes.inverted()

        start, end, si, ei, chr, orientation, span = ext
        flank = span / scale / 2
        xstart, xend = x - flank, x + flank

        cv = lambda t: xstart + abs(t - startbp) / scale

        # Chromosome
        ax.plot((xstart, xend), (y, y), color="gray", transform=tr, \
                lw=2, zorder=1)

        genes = bed[si: ei + 1]
        startbp, endbp = start.start, end.end
        if orientation == '-':
            startbp, endbp = endbp, startbp

        if switch:
            chr = switch.get(chr, chr)
        label = "-".join((human_size(startbp)[:-2], human_size(endbp)))

        height = .012
        self.gg = {}
        # Genes
        for g in genes:
            gstart, gend = g.start, g.end
            strand = g.strand
            if strand == '-':
                gstart, gend = gend, gstart
            if orientation == '-':
                strand = "+" if strand == "-" else "-"

            x1, x2 = cv(gstart), cv(gend)
            a, b = tr.transform((x1, y)), tr.transform((x2, y))
            a, b = inv.transform(a), inv.transform(b)
            self.gg[g.accn] = (a, b)
            #ax.text(a[0], a[1], g.accn[-6:], size=6, rotation=lr + 90)

            color = "b" if strand == "+" else "g"
            gp = Glyph(ax, x1, x2, y, height, gradient=False, fc=color, zorder=2)
            gp.set_transform(tr)

        ha, va = layout.ha, layout.va
        if ha == "left":
            xx = xstart
            ha = "right"
        elif ha == "right":
            xx = xend
            ha = "left"
        else:
            xx = x
            ha = "center"

        if va == "top":
            yy = y + pad
        elif va == "bottom":
            yy = y - pad
        else:
            yy = y

        l = np.array((xx, yy))
        trans_angle = ax.transAxes.transform_angles(np.array((lr, )),
                                                    l.reshape((1, 2)))[0]
        p3 = pad / 3
        lx, ly = l
        ax.text(lx, ly, chr + "\n ", color=layout.color,
                    ha=ha, va="center", size=10, rotation=trans_angle)
        ax.text(lx, ly, _(" \n \n" + label), color="k",
                    ha=ha, va="center", size=10, rotation=trans_angle)



def main():
    p = OptionParser(__doc__)
    p.add_option("--switch",
                 help="Rename the seqid with two-column file [default: %default]")
    opts, args, iopts = set_image_options(p, figsize="6x6")

    if len(args) != 3:
        sys.exit(not p.print_help())

    datafile, bedfile, layoutfile = args
    bed = Bed(bedfile)
    order = bed.order
    bf = BlockFile(datafile)
    lo = Layout(layoutfile)
    switch = opts.switch
    switch = DictFile(switch, delimiter="\t") if switch else None

    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    exts = []
    for i in xrange(bf.ncols):
        ext = bf.get_extent(i, order)
        exts.append(ext)

    maxspan = max(exts, key=lambda x: x[-1])[-1]
    scale = maxspan / .65

    gg = {}
    ymids = []
    for i in xrange(bf.ncols):
        ext = exts[i]
        r = Region(root, ext, lo[i], bed, scale, switch)
        gg.update(r.gg)
        ymids.append(r.y)

    for i, j in lo.edges:
        for ga, gb in bf.iter_pairs(i, j):
            a, b = gg[ga], gg[gb]
            ymid = (ymids[i] + ymids[j]) / 2
            Shade(root, a, b, ymid)

        for ga, gb in bf.iter_pairs(i, j, highlight=True):
            a, b = gg[ga], gg[gb]
            ymid = (ymids[i] + ymids[j]) / 2
            Shade(root, a, b, ymid, highlight=True)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


if __name__ == '__main__':
    main()
