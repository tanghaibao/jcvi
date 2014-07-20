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

from jcvi.compara.synteny import BlockFile
from jcvi.formats.bed import Bed
from jcvi.formats.base import DictFile
from jcvi.utils.cbook import human_size
from jcvi.apps.base import OptionParser

from jcvi.graphics.glyph import Glyph, RoundLabel
from jcvi.graphics.base import mpl, plt, savefig, markup, \
            Path, PathPatch, AbstractLayout


class LayoutLine (object):

    def __init__(self, row, delimiter=','):
        self.hidden = row[0] == '*'
        if self.hidden:
            row = row[1:]
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.x = float(args[0])
        self.y = float(args[1])
        self.rotation = int(args[2])
        self.ha = args[3]
        self.va = args[4]
        self.color = args[5]
        self.ratio = 1
        if len(args) > 6:
            self.ratio = float(args[6])


class Layout (AbstractLayout):

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

        self.assign_colors()


class Shade (object):

    def __init__(self, ax, a, b, ymid, highlight=False, ec="k", fc="k",
                    alpha=.2, lw=1, zorder=1):
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
            ec = fc = highlight

        pp = PathPatch(path, ec=ec, fc=fc, alpha=alpha,
                     lw=lw, zorder=zorder)
        ax.add_patch(pp)


class Region (object):

    def __init__(self, ax, ext, layout, bed, scale, switch=None, chr_label=True,
                 pad=.04, vpad=.012):
        x, y = layout.x, layout.y
        ratio = layout.ratio
        scale /= ratio
        self.y = y
        lr = layout.rotation
        tr = mpl.transforms.Affine2D().\
                    rotate_deg_around(x, y, lr) + ax.transAxes
        inv = ax.transAxes.inverted()

        start, end, si, ei, chr, orientation, span = ext
        flank = span / scale / 2
        xstart, xend = x - flank, x + flank
        self.xstart, self.xend = xstart, xend

        cv = lambda t: xstart + abs(t - startbp) / scale
        hidden = layout.hidden

        # Chromosome
        if not hidden:
            ax.plot((xstart, xend), (y, y), color="gray", transform=tr, \
                    lw=2, zorder=1)

        self.genes = genes = bed[si: ei + 1]
        startbp, endbp = start.start, end.end
        if orientation == '-':
            startbp, endbp = endbp, startbp

        if switch:
            chr = switch.get(chr, chr)
        label = "-".join((human_size(startbp, target="Mb")[:-2],
                          human_size(endbp, target="Mb")))

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

            color = "b" if strand == "+" else "g"
            if not hidden:
                gp = Glyph(ax, x1, x2, y, height, gradient=False, fc=color, zorder=3)
                gp.set_transform(tr)

        ha, va = layout.ha, layout.va

        hpad = .02
        if ha == "left":
            xx = xstart - hpad
            ha = "right"
        elif ha == "right":
            xx = xend + hpad
            ha = "left"
        else:
            xx = x
            ha = "center"

        # Tentative solution to labels stick into glyph
        magic = 40.
        cc = abs(lr) / magic if abs(lr) > magic else 1
        if va == "top":
            yy = y + cc * pad
        elif va == "bottom":
            yy = y - cc * pad
        else:
            yy = y

        l = np.array((xx, yy))
        trans_angle = ax.transAxes.transform_angles(np.array((lr, )),
                                                    l.reshape((1, 2)))[0]
        lx, ly = l
        if not hidden and chr_label:
            ax.text(lx, ly + vpad, markup(chr), color=layout.color,
                        ha=ha, va="center", rotation=trans_angle)
            ax.text(lx, ly - vpad, label, color="k",
                        ha=ha, va="center", rotation=trans_angle)


class Synteny (object):

    def __init__(self, fig, root, datafile, bedfile, layoutfile,
                 switch=None, tree=None, chr_label=True, pad=.04):

        w, h = fig.get_figwidth(), fig.get_figheight()
        bed = Bed(bedfile)
        order = bed.order
        bf = BlockFile(datafile)
        lo = Layout(layoutfile)
        switch = DictFile(switch, delimiter="\t") if switch else None

        exts = []
        for i in xrange(bf.ncols):
            ext = bf.get_extent(i, order)
            exts.append(ext)

        maxspan = max(exts, key=lambda x: x[-1])[-1]
        scale = maxspan / .65

        self.gg = gg = {}
        self.rr = []
        ymids = []
        vpad = .012 * w / h
        for i in xrange(bf.ncols):
            ext = exts[i]
            r = Region(root, ext, lo[i], bed, scale, switch, chr_label=chr_label, vpad=vpad)
            self.rr.append(r)
            # Use tid and accn to store gene positions
            gg.update(dict(((i, k), v) for k, v in r.gg.items()))
            ymids.append(r.y)

        for i, j in lo.edges:
            for ga, gb, h in bf.iter_pairs(i, j):
                a, b = gg[(i, ga)], gg[(j, gb)]
                ymid = (ymids[i] + ymids[j]) / 2
                Shade(root, a, b, ymid, fc="gainsboro", lw=0, alpha=1)

            for ga, gb, h in bf.iter_pairs(i, j, highlight=True):
                a, b = gg[(i, ga)], gg[(j, gb)]
                ymid = (ymids[i] + ymids[j]) / 2
                Shade(root, a, b, ymid, alpha=1, highlight=h, zorder=2)

        if tree:
            from jcvi.graphics.tree import draw_tree, read_trees

            trees = read_trees(tree)
            ntrees = len(trees)
            logging.debug("A total of {0} trees imported.".format(ntrees))
            xiv = 1. / ntrees
            yiv = .3
            xstart = 0
            ystart = min(ymids) - .4
            for i in xrange(ntrees):
                ax = fig.add_axes([xstart, ystart, xiv, yiv])
                label, outgroup, tx = trees[i]
                draw_tree(ax, tx, outgroup=outgroup, rmargin=.4)
                xstart += xiv
                RoundLabel(ax, .5, .3, label, fill=True, fc="lavender", color="r")


def draw_gene_legend(ax, x1, x2, ytop):
    d = .04
    ax.plot([x1, x1 + d], [ytop, ytop], "b:", lw=2)
    ax.plot([x1 + d], [ytop], "b>", mec="b")
    ax.plot([x2, x2 + d], [ytop, ytop], "g:", lw=2)
    ax.plot([x2], [ytop], "g<", mec="g")


def main():
    p = OptionParser(__doc__)
    p.add_option("--switch",
                 help="Rename the seqid with two-column file [default: %default]")
    p.add_option("--tree",
                 help="Display trees on the bottom of the figure [default: %default]")
    opts, args, iopts = p.set_image_options(figsize="8x7")

    if len(args) != 3:
        sys.exit(not p.print_help())

    datafile, bedfile, layoutfile = args
    switch = opts.switch
    tree = opts.tree

    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Synteny(fig, root, datafile, bedfile, layoutfile,
            switch=switch, tree=tree)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
