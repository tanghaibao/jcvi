#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Gradient gene features
"""

import os.path as op
import sys

import numpy as np
from random import choice, shuffle, random, randint

from jcvi.apps.base import OptionParser, ActionDispatcher
from jcvi.graphics.base import plt, Rectangle, CirclePolygon, Polygon, \
            savefig, get_map


tstep = .05
Timing = np.arange(0, 1 + tstep, tstep)
arrowprops = dict(arrowstyle="fancy", fc="lightslategray", ec="lightslategray",
            connectionstyle="arc3,rad=-0.05")


class Bezier (object):
    """
    Cubic bezier curve, see the math:
    <http://www.moshplant.com/direct-or/bezier/math.html>
    p0 : origin, p1, p2 :control, p3: destination
    """
    def __init__(self, ax, p0, p1, p2, p3, color='m', alpha=.2):
        pts = (p0, p1, p2, p3)
        px, py = zip(*pts)
        xt = self.get_array(px)
        yt = self.get_array(py)

        ax.plot(xt, yt, "-", color=color, alpha=alpha)

    def get_array(self, pts, t=Timing):
        p0, p1, p2, p3 = pts

        # Get the coeffiencients
        c = 3 * (p1 - p0)
        b = 3 * (p2 - p1) - c
        a = p3 - p0 - c - b

        tsquared = t ** 2
        tcubic = tsquared * t
        return a * tcubic + b * tsquared + c * t + p0


class RoundLabel (object):
    """Round rectangle around the text label
    """
    def __init__(self, ax, x1, x2, t, lw=0, fill=False, fc="lavender", **kwargs):

        ax.text(x1, x2, t, ha="center",
            bbox=dict(boxstyle="round", fill=fill, fc=fc, lw=lw), **kwargs)


class RoundRect (object):
    """Round rectangle directly
    """
    def __init__(self, ax, xy, width, height, shrink=.1, label=None, **kwargs):

        shrink *= height
        x, y= xy
        pts = []
        # plot the four rounded cap one by one
        pts += plot_cap((x + width - shrink, y + height - shrink),
                             np.radians(range(0, 90)), shrink)
        pts += [[x + width - shrink, y + height], [x + shrink, y + height]]
        pts += plot_cap((x + shrink, y + height - shrink),
                             np.radians(range(90, 180)), shrink)
        pts += [[x, y + height - shrink], [x, y + shrink]]
        pts += plot_cap((x + shrink, y + shrink),
                             np.radians(range(180, 270)), shrink)
        pts += [[x + shrink, y], [x + width - shrink, y]]
        pts += plot_cap((x + width - shrink, y + shrink),
                             np.radians(range(270, 360)), shrink)
        pts += [[x + width, y + shrink], [x + width, y + height - shrink]]
        p1 = Polygon(pts, **kwargs)
        ax.add_patch(p1)
        # add a white transparency ellipse filter
        if label:
            ax.text(x + width / 2,y + height / 2,label, size=10,
                      ha="center", va="center", color="w")


class DoubleSquare (object):
    """Square with a double-line margin
    """
    def __init__(self, ax, x, y, radius=.01, **kwargs):

        d = radius * 1.5
        ax.add_patch(Rectangle((x - d, y - d), 2 * d, 2 * d,
                     fc="w", ec="k", zorder=10))
        d = radius
        ax.add_patch(Rectangle((x - d, y - d), 2 * d, 2 * d,
                     zorder=10, **kwargs))


class DoubleCircle (object):
    """Circle with a double-line margin
    """
    def __init__(self, ax, x, y, radius=.01, **kwargs):

        ax.add_patch(CirclePolygon((x, y), radius * 1.4,
                     resolution=50, fc="w", ec="k"))
        ax.add_patch(CirclePolygon((x, y), radius,
                     resolution=50, **kwargs))


class TextCircle (object):
    """Circle with a character wrapped in
    """
    def __init__(self, ax, x, y, label, radius=.02, fc="k", color="w",
                size=12, zorder=4, **kwargs):

        circle = CirclePolygon((x, y), radius, resolution=20, \
                                fc=fc, ec=fc, zorder=zorder, **kwargs)
        ax.add_patch(circle)
        ax.text(x, y, label, ha="center", va="center", color=color,
                size=size, zorder=zorder + 1, **kwargs)


class BaseGlyph (list):

    def __init__(self, ax):
        self.ax = ax

    def add_patches(self):
        for p in self:
            self.ax.add_patch(p)

    def set_transform(self, tr):
        for p in self:
            p.set_transform(tr)


class Glyph (BaseGlyph):
    """Draws gradient rectangle
    """
    def __init__(self, ax, x1, x2, y, height=.03, gradient=True, fc="gray", **kwargs):

        super(Glyph, self).__init__(ax)
        width = x2 - x1
        # Frame around the gradient rectangle
        p1 = (x1, y - .5 * height)
        self.append(Rectangle(p1, width, height, fc=fc,
            lw=0, **kwargs))

        # Several overlaying patches
        if gradient:
            for cascade in np.arange(.1, .55, .05):
                p1 = (x1, y - height * cascade)
                self.append(Rectangle(p1, width, 2 * cascade * height,
                    fc='w', lw=0, alpha=.1))

        self.add_patches()


class ExonGlyph (BaseGlyph):
    """Multiple rectangles linked together.
    """
    def __init__(self, ax, x, y, mrnabed, exonbeds, height=.03, ratio=1,
                 align="left", **kwargs):

        super(ExonGlyph, self).__init__(ax)
        start, end = mrnabed.start, mrnabed.end
        xa = lambda a: x + (a - start) * ratio
        xb = lambda a: x - (end - a) * ratio
        xc = xa if align == "left" else xb

        Glyph(ax, xc(start), xc(end), y, height=height / 3)
        for b in exonbeds:
            bstart, bend = b.start, b.end
            Glyph(ax, xc(bstart), xc(bend), y, fc="orange")


class GeneGlyph (BaseGlyph):
    """Draws an oriented gene symbol, with color gradient, to represent genes
    """
    def __init__(self, ax, x1, x2, y, height, gradient=True, tip=.0025, \
                 color="k", shadow=False, **kwargs):

        super(GeneGlyph, self).__init__(ax)
        # Figure out the polygon vertices first
        orientation = 1 if x1 < x2 else -1
        level = 10
        tip = min(tip, abs(x1 - x2))
        # Frame
        p1 = (x1, y - height * .5)
        p2 = (x2 - orientation * tip, y - height * .5)
        p3 = (x2, y)
        p4 = (x2 - orientation * tip, y + height * .5)
        p5 = (x1, y + .5*height)
        if "fc" not in kwargs:
            kwargs["fc"] = color
        if "ec" not in kwargs:
            kwargs["ec"] = color
        P = Polygon([p1, p2, p3, p4, p5], **kwargs)
        self.append(P)

        if gradient:
            zz = kwargs.get("zorder", 1)
            zz += 1
            # Patch (apply white mask)
            for cascade in np.arange(0, .5, .5 / level):
                p1 = (x1, y - height * cascade)
                p2 = (x2 - orientation * tip, y - height * cascade)
                p3 = (x2, y)
                p4 = (x2 - orientation * tip, y + height * cascade)
                p5 = (x1, y + height * cascade)
                self.append(Polygon([p1, p2, p3, p4, p5], fc='w', \
                        lw=0, alpha=.2, zorder=zz))

        if shadow:
            import matplotlib.patheffects as pe
            P.set_path_effects([pe.withSimplePatchShadow((1, -1), \
                                alpha=.4, patch_alpha=1)])

        self.add_patches()


class CartoonRegion (object):
    """
    Draw a collection of GeneGlyphs along chromosome.
    """
    def __init__(self, n, k=12):
        # Chromosome
        self.n = n
        self.orientations = [choice([-1, 1]) for i in xrange(n)]
        self.colors = self.assign_colors(k)

    def draw(self, ax, x, y, gene_len=.012, strip=True, color=True):
        if strip:
            self.strip()

        t = gene_len * 1.2
        length = t * (self.n + 1)
        x1, x2 = x - length / 2, x + length / 2
        self.x1, self.x2 = x1, x2
        self.y = y
        ax.plot((x1, x2), (y, y), color="gray", lw=2, zorder=1)
        bit = .008
        xs = (x1 - 2 * bit, x1 - bit, x2 + bit, x2 + 2 * bit)
        ax.plot(xs, [y] * 4, ".", lw=2, color="gray")
        pos = np.arange(x1 + t, x2, t)[:self.n]
        assert len(pos) == self.n, "len(pos) = {0}".format(len(pos))

        gl = gene_len / 2
        for x, c, o in zip(pos, self.colors, self.orientations):
            x1, x2 = x - gl, x + gl
            if o < 0:
                x1, x2 = x2, x1
            if not color and c != 'k':
                c = 'w'
            GeneGlyph(ax, x1, x2, y, gene_len, color=c, ec='k',
                      gradient=False, shadow=True, zorder=10)

    def assign_colors(self, k):
        from matplotlib.colors import rgb2hex

        colorset = get_map('Paired', 'qualitative', k).mpl_colors
        colorset = [rgb2hex(x) for x in colorset]
        cs = colorset + ['w'] * (self.n - k - 1)
        shuffle(cs)
        return cs[:self.n / 2] + ['k'] + cs[self.n / 2:]

    def delete(self, p, waiver=None):
        if waiver and self.colors[p] in waiver:
            return
        self.colors.pop(p)
        self.orientations.pop(p)
        self.n -= 1

    def insert(self, p):
        self.colors.insert(p, 'w')
        self.orientations.insert(p, choice([-1, 1]))
        self.n += 1

    def truncate(self, b, e):
        b = max(b, 0)
        e = min(self.n, e)
        self.colors = self.colors[b:e]
        self.orientations = self.orientations[b:e]
        self.n = e - b

    def assign_flankers(self):
        lf, p, rf = self.find_k()
        self.flanks = [self.colors[lf], self.colors[rf]]
        return p

    def truncate_between_flankers(self, target=0):
        try:
            lf, rf = self.flanks
        except:
            self.assign_flankers()
            lf, rf = self.flanks
        lf = self.colors.index(lf) if lf in self.colors else -1
        rf = self.colors.index(rf) if rf in self.colors else -1
        assert lf >= 0 or rf >= 0
        if rf < 0:
            rf = lf
        if lf < 0:
            lf = rf
        if rf + 1 - lf < target:
            gap = target - rf - 1 + lf
            lf -= gap / 2
            rf += gap / 2
        self.truncate(lf, rf + 1)

    def strip(self):
        while self.colors[0] == 'w':
            self.delete(0)
        while self.colors[-1] == 'w':
            self.delete(self.n - 1)

    def find_k(self):
        p = self.colors.index('k')
        lf = max(i for i, c in enumerate(self.colors[:p]) if c != 'w')
        rf = min(i for i, c in enumerate(self.colors[p + 1:]) if c != 'w')
        return lf, p, rf + p + 1

    def evolve(self, mode='S', target=10):
        n = self.n
        assert mode in ('S', 'F', 'G')
        keep_k = mode == 'S'
        p = self.assign_flankers()
        waiver = self.flanks[:]
        if mode == 'S':
            waiver += ['k']
        if mode == 'F':
            self.delete(p)
        elif mode == 'G':
            left_score = sum(1 for x in self.colors[:p] if x != 'w')
            right_score = sum(1 for x in self.colors[p + 1:] if x != 'w')
            if left_score > right_score:
                self.colors[:p] = ['w'] * p
            else:
                self.colors[p + 1:] = ['w'] * (self.n - p - 1)
        while self.nonwhites > target:
            if random() > .35:
                self.delete(randint(0, self.n - 1), waiver=waiver)
            if random() > .65 and self.n < n * .8:
                self.insert(randint(0, self.n - 1))

    @property
    def nonwhites(self):
        return sum(1 for x in self.colors if x != 'w')


def plot_cap(center, t, r):
    x, y = center
    return zip(x + r * np.cos(t), y + r * np.sin(t))


def main():

    actions = (
        ('demo', 'run a demo to showcase some common usages of various glyphs'),
        ('gff', 'draw exons for genes based on gff files'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_cds_beds(gffile, noUTR=False):
    from jcvi.formats.gff import Gff

    mrnabed = None
    cdsbeds = []
    gf = Gff(gffile)
    for g in gf:
        if g.type == "mRNA":
            mrnabed = g.bedline
        elif g.type == "CDS":
            cdsbeds.append(g.bedline)

    if noUTR:
        mrnabed.start = min(x.start for x in cdsbeds)
        mrnabed.end = max(x.end for x in cdsbeds)

    return mrnabed, cdsbeds


def get_setups(gffiles, canvas=.6, noUTR=False):
    setups = []
    for gffile in gffiles:
        genename = op.basename(gffile).rsplit(".", 1)[0]
        mrnabed, cdsbeds = get_cds_beds(gffile, noUTR=noUTR)
        setups.append((genename, mrnabed, cdsbeds))

    genenames, mrnabeds, cdsbedss = zip(*setups)
    maxspan = max(x.span for x in mrnabeds)
    ratio = canvas / maxspan
    return setups, ratio


def gff(args):
    """
    %prog gff *.gff

    Draw exons for genes based on gff files. Each gff file should contain only
    one gene, and only the "mRNA" and "CDS" feature will be drawn on the canvas.
    """
    align_choices = ("left", "center", "right")
    p = OptionParser(gff.__doc__)
    p.add_option("--align", default="left", choices=align_choices,
                 help="Horizontal alignment [default: %default]")
    p.add_option("--noUTR", default=False, action="store_true",
                 help="Do not plot UTRs [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fig = plt.figure(1, (8, 5))
    root = fig.add_axes([0, 0, 1, 1])

    gffiles = args
    ngenes = len(gffiles)

    canvas = .6
    setups, ratio = get_setups(gffiles, canvas=canvas, noUTR=opts.noUTR)
    align = opts.align
    xs = .2 if align == "left" else .8
    yinterval = canvas / ngenes
    ys = .8
    tip = .01
    for genename, mrnabed, cdsbeds in setups:
        ExonGlyph(root, xs, ys, mrnabed, cdsbeds, ratio=ratio, align=align)
        if align == "left":
            root.text(xs - tip, ys, genename, ha="right", va="center")
        elif align == "right":
            root.text(xs + tip, ys, genename, ha="left", va="center")
        ys -= yinterval

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = "exons.pdf"
    savefig(figname, dpi=300)


def demo(args):
    """
    %prog demo

    Draw sample gene features to illustrate the various fates of duplicate
    genes - to be used in a book chapter.
    """
    p = OptionParser(demo.__doc__)
    opts, args = p.parse_args(args)

    fig = plt.figure(1, (8, 5))
    root = fig.add_axes([0, 0, 1, 1])

    panel_space = .23
    dup_space = .025
    # Draw a gene and two regulatory elements at these arbitrary locations
    locs = [(.5, .9), # ancestral gene
            (.5, .9 - panel_space + dup_space), # identical copies
            (.5, .9 - panel_space - dup_space),
            (.5, .9 - 2 * panel_space + dup_space), # degenerate copies
            (.5, .9 - 2 * panel_space - dup_space),
            (.2, .9 - 3 * panel_space + dup_space), # sub-functionalization
            (.2, .9 - 3 * panel_space - dup_space),
            (.5, .9 - 3 * panel_space + dup_space), # neo-functionalization
            (.5, .9 - 3 * panel_space - dup_space),
            (.8, .9 - 3 * panel_space + dup_space), # non-functionalization
            (.8, .9 - 3 * panel_space - dup_space),
            ]

    default_regulator = "gm"
    regulators = [default_regulator,
            default_regulator, default_regulator,
            "wm", default_regulator,
            "wm", "gw",
            "wb", default_regulator,
            "ww", default_regulator,
            ]

    width = .24
    for i, (xx, yy) in enumerate(locs):
        regulator = regulators[i]
        x1, x2 = xx - .5 * width, xx + .5 * width
        Glyph(root, x1, x2, yy)
        if i == 9:  # upper copy for non-functionalization
            continue

        # coding region
        x1, x2 = xx - .16 * width, xx + .45 * width
        Glyph(root, x1, x2, yy, fc="k")

        # two regulatory elements
        x1, x2 = xx - .4 * width, xx - .28 * width
        for xx, fc in zip((x1, x2), regulator):
            if fc == 'w':
                continue

            DoubleCircle(root, xx, yy, fc=fc)

        rotation = 30
        tip = .02
        if i == 0:
            ya = yy + tip
            root.text(x1, ya, "Flower", rotation=rotation, va="bottom")
            root.text(x2, ya, "Root", rotation=rotation, va="bottom")
        elif i == 7:
            ya = yy + tip
            root.text(x2, ya, "Leaf", rotation=rotation, va="bottom")

    # Draw arrows between panels (center)
    arrow_dist = .08
    ar_xpos = .5
    for ar_ypos in (.3, .53, .76):
        root.annotate(" ", (ar_xpos, ar_ypos),
                (ar_xpos, ar_ypos + arrow_dist),
                arrowprops=arrowprops)

    ar_ypos = .3
    for ar_xpos in (.2, .8):
        root.annotate(" ", (ar_xpos, ar_ypos),
                (.5, ar_ypos + arrow_dist),
                arrowprops=arrowprops)

    # Duplication, Degeneration
    xx = .6
    ys = (.76, .53)
    processes = ("Duplication", "Degeneration")
    for yy, process in zip(ys, processes):
        root.text(xx, yy + .02, process, fontweight="bold")

    # Label of fates
    xs = (.2, .5, .8)
    fates = ("Subfunctionalization", "Neofunctionalization",
            "Nonfunctionalization")
    yy = .05
    for xx, fate in zip(xs, fates):
        RoundLabel(root, xx, yy, fate)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = "demo.pdf"
    savefig(figname, dpi=300)


if __name__ == '__main__':
    main()
