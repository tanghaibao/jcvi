#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog demo

Illustrate three different types of alignments.
- Pairwise sequence alignment, aka, "dot plot"
- Read alignment, similar to the visualization of a BAM file
- Optical map alignment, matchings between restriction fragments
"""


import sys

from random import randint, choice

from jcvi.utils.iter import pairwise
from jcvi.utils.range import range_overlap
from jcvi.graphics.chromosome import Chromosome, HorizontalChromosome
from jcvi.graphics.glyph import GeneGlyph
from jcvi.graphics.base import Rectangle, plt, savefig, normalize_axes, set2
from jcvi.apps.base import OptionParser


class BaseAlign (object):

    def __init__(self, fig, xywh, xpad=0, ypad=0, xmax=100):
        x, y, w, h = xywh
        self.ax = fig.add_axes(xywh)
        self.sax = fig.add_axes([x + xpad * w, y + ypad * h,
                                (1 - 2 * xpad) * w, (1 - 2 * ypad) * h])
        self.amax = self.bmax = xmax
        self.a = [(1, xmax)]
        self.b = [(1, xmax)]
        self.apatch = self.bpatch = None
        self.xpad = xpad
        self.ypad = ypad
        self.canvas = 1 - 2 * xpad

    def convert(self, pos, xmax):
        return self.xpad + pos * self.canvas / xmax

    def invert(self, a, b):
        self.a = [(1, a), (a, b), (b, self.amax)]
        self.b = [(1, a), (b, a), (b, self.bmax)]
        self.apatch = (self.convert(a, self.amax),
                       self.convert(b, self.amax))
        self.bpatch = (self.convert(a, self.bmax),
                       self.convert(b, self.bmax))

    def delete(self, a, b):
        self.bmax -= b - a
        self.a = [(1, a), (b, self.amax)]
        self.b = [(1, a), (a, self.bmax)]
        self.apatch = (self.convert(a, self.amax),
                       self.convert(b, self.amax))

    def duplicate(self, a, b, gap=0):
        self.bmax += b - a + gap
        self.a = [(1, b), (a, self.amax)]
        self.b = [(1, b), (b + gap, self.bmax)]
        self.apatch = (self.convert(a, self.amax),
                       self.convert(b, self.amax))
        self.bpatch = (self.convert(a, self.bmax),
                       self.convert(b, self.bmax),
                       self.convert(b + gap, self.bmax),
                       self.convert(2 * b - a + gap, self.bmax))


class PairwiseAlign (BaseAlign):

    def __init__(self, fig, xywh, xpad=.15, ypad=.15):
        super(PairwiseAlign, self).__init__(fig, xywh, xpad, ypad)

    def draw(self, width=.03):
        HorizontalChromosome(self.ax, self.xpad, 1 - self.xpad,
                             self.ypad - .05, height=width * 1.5,
                             patch=self.apatch, lw=2)
        Chromosome(self.ax, self.xpad - .05, self.ypad, 1 - self.ypad,
                   width=width, patch=self.bpatch, lw=2)
        for a, b in zip(self.a, self.b):
            self.sax.plot(a, b, "-", color="darkslategrey", lw=2)
        self.sax.set_xticklabels([])
        self.sax.set_yticklabels([])
        self.sax.set_xlim((1, self.amax))
        self.sax.set_ylim((1, self.bmax))
        normalize_axes(self.ax)


class ReadAlign (BaseAlign):

    def __init__(self, fig, xywh, xpad=.05, ypad=.2, readlen=6, gap=3):
        super(ReadAlign, self).__init__(fig, xywh, xpad, ypad)
        self.readlen = readlen
        self.gap = gap
        self.reads = []
        self.ymax = 12
        self.ntracks = 0
        self.layout(1, self.amax)

    def layout(self, start, end, maxtracks=8):
        readrange = 2 * self.readlen + self.gap
        end -= readrange
        assert start < end, "end must be > start + readlen"
        reads = []
        for x in xrange(100):
            pos = randint(start, end)
            reads.append(PairedRead(pos, readlen=self.readlen, gap=self.gap))
        reads, ntracks = self.arrange(reads, self.ntracks, maxtracks=maxtracks)
        self.reads += reads
        self.ntracks += ntracks

    def arrange(self, reads, ntracks, maxtracks=8):
        track_ends = [0]
        reads.sort(key=lambda x: x.start)
        for r in reads:
            m = min(track_ends)
            mi = track_ends.index(m)
            if r.start > m + .005:
                track_ends[mi] = r.end
            else:
                if len(track_ends) >= maxtracks:
                    continue
                track_ends.append(r.end)
                mi = len(track_ends) - 1
            r.set_y(ntracks + mi)
        ntracks = len(track_ends)
        reads = [x for x in reads if x.y is not None]
        return reads, ntracks

    def remove(self, a, b, maxtracks=0):
        self.reads = [r for r in self.reads \
                      if not (a <= r.start <= b and a <= r.end <= b \
                              and r.y >= maxtracks)]

    def draw(self, width=.03):
        HorizontalChromosome(self.ax, self.xpad, 1 - self.xpad,
                             self.ypad - width / 2, height=width * 1.5,
                             patch=self.apatch, lw=2)
        for r in self.reads:
            r.draw(self.sax)
        self.sax.set_xlim((1, self.amax))
        self.sax.set_ylim((-1, self.ymax))
        normalize_axes(self.ax)
        self.sax.set_axis_off()

    def highlight(self, a, b):
        self.apatch = (self.convert(a, self.amax),
                       self.convert(b, self.amax))
        self.sax.plot((a, a), (-1, self.ntracks), "m-", lw=2)
        self.sax.plot((b, b), (-1, self.ntracks), "m-", lw=2)

    def invert(self, a, b):
        reads = []
        for r in self.reads:
            r.set_y(None)
            keep = True
            if r.start < a < r.end or r.start < b < r.end:
                adist, bdist = abs(a - r.mid), abs(b - r.mid)
                flipr = r.r2 if adist > bdist else r.r1
                flipr.x1 = a + b - flipr.x1
                flipr.x2 = a + b - flipr.x2
                flipr.color = 'y'
                if adist > self.gap and bdist > self.gap:
                    keep = False
            if keep:
                reads.append(r)
        self.reads, self.ntracks = self.arrange(reads, 0)
        self.highlight(a, b)

    def delete(self, a, b):
        self.remove(a, b)
        for r in self.reads:
            r.breakpoint(a, 'g', 'lightgrey')
            r.breakpoint(b, 'lightgrey', 'g')
        self.highlight(a, b)

    def duplicate(self, a, b):
        self.layout(1, self.amax, maxtracks=4)
        self.remove(1, a, maxtracks=6)
        self.remove(b, self.amax, maxtracks=6)
        for r in self.reads:
            r.paint(a, b, 'r')
            r.breakpoint(a, 'k', 'r')
            r.breakpoint(b, 'r', 'k')
            r.breakpoint(a, 'lightgrey', 'r', ystart=6)
            r.breakpoint(b, 'r', 'lightgrey', ystart=6)
        self.highlight(a, b)


class OpticalMapAlign (BaseAlign):

    def __init__(self, fig, xywh, xpad=.05, ypad=.3):
        super(OpticalMapAlign, self).__init__(fig, xywh, xpad, ypad)
        self.om1 = self.from_silico()
        self.om2 = self.om1[:]

    def make_wiggles(self, length, wiggle=3):
        ar = [0]
        while len(ar) <= length:
            nm = ar[-1] + choice((-1, 1))
            if not 0 <= nm < wiggle:
                continue
            ar.append(nm)
        return ar

    def make_colors(self, length, palette=set2[:4]):
        ar = []
        while len(ar) <= length:
            nm = choice(palette)
            if ar and nm == ar[-1]:
                continue
            ar.append(nm)
        return ar

    def from_silico(self, filename="Ecoli.silico", nfrags=40):
        fp = open(filename)
        fp.next()
        frags = [0] + [int(x) for x in fp.next().split()]
        return frags[:nfrags]

    def draw_track(self, ar, ystart=0):
        nfrags = len(ar)
        wiggles = self.make_wiggles(nfrags)
        colors = self.make_colors(nfrags)
        j = 0
        for a, b in pairwise(ar):
            if b - a < max(ar) / 200:
                continue
            j += 1
            wiggle, color = wiggles[j], colors[j]
            yf = ystart + wiggle * 1. / max(wiggles)
            p = Rectangle((a, yf), b - a, 1, color=color)
            self.sax.add_patch(p)

    def draw(self):
        self.draw_track(self.om1)
        self.draw_track(self.om2, ystart=10)
        self.sax.set_xlim((-1, max(self.om1)))
        self.sax.set_ylim((-1, 12))
        normalize_axes(self.ax)
        self.sax.set_axis_off()


class SingleRead (object):

    def __init__(self, start, readlen, sign=1):
        self.x1 = start
        self.x2 = start + sign * readlen
        self.y = None
        self.color = 'k'
        self.broken = None

    @property
    def sign(self):
        return 1 if self.x2 >= self.x1 else -1

    @property
    def start(self):
        return min(self.x1, self.x2)

    @property
    def end(self):
        return max(self.x1, self.x2)

    @property
    def span(self):
        return self.end - self.start + 1

    def draw(self, ax, height=.6):
        if self.broken is None:
            GeneGlyph(ax, self.x1, self.x2, self.y, height, tip=2,
                      color=self.color, gradient=True)
        else:
            a, lcolor, rcolor = self.broken
            if self.sign < 0:
                lcolor, rcolor = rcolor, lcolor
            GeneGlyph(ax, self.x1, a, self.y, height, tip=0,
                      color=lcolor, gradient=True)
            GeneGlyph(ax, a, self.x2, self.y, height, tip=2,
                      color=rcolor, gradient=True)

    def breakpoint(self, a, lcolor, rcolor):
        if a > self.end:
            self.color = lcolor
        elif a < self.start:
            self.color = rcolor
        else:
            self.broken = (a, lcolor, rcolor)


class PairedRead (object):

    def __init__(self, start, readlen, gap):
        self.r1 = SingleRead(start, readlen)
        self.r2 = SingleRead(start + gap + 2 * readlen, readlen, sign=-1)
        self.color = 'k'
        self.y = None

    @property
    def start(self):
        return min(self.r1.start, self.r2.start)

    @property
    def end(self):
        return max(self.r1.end, self.r2.end)

    @property
    def i1(self):
        return min(self.r1.end, self.r2.end)

    @property
    def i2(self):
        return max(self.r1.start, self.r2.start)

    @property
    def mid(self):
        return (self.start + self.end) * .5

    def set_y(self, y):
        self.y = y
        self.r1.y = self.r2.y = y

    def draw(self, ax):
        self.r1.draw(ax)
        self.r2.draw(ax)
        ax.plot((self.i1, self.i2), (self.y, self.y), "-",
                 color=self.color)

    def paint(self, a, b, color):
        if range_overlap((0, self.start + 1 , self.end - 1),
                         (0, a, b)):
            self.r1.color = self.r2.color = self.color = color

    def breakpoint(self, a, lcolor, rcolor, ystart=0):
        if not self.start < a < self.end:
            return
        if self.y < ystart:
            return
        self.color = lcolor if a > self.mid else rcolor
        self.r1.breakpoint(a, lcolor, rcolor)
        self.r2.breakpoint(a, lcolor, rcolor)


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = p.set_image_options(figsize="9x7")

    if len(args) != 1:
        sys.exit(not p.print_help())

    mode, = args
    assert mode == "demo"

    w = .33
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    p = PairwiseAlign(fig, [0, 2 * w, w, w])
    p.invert(30, 70)
    p.draw()

    p = PairwiseAlign(fig, [0, w, w, w])
    p.delete(30, 70)
    p.draw()

    p = PairwiseAlign(fig, [0, 0, w, w])
    p.duplicate(30, 70, gap=5)
    p.draw()

    p = ReadAlign(fig, [w, 2 * w, w, w])
    p.invert(30, 70)
    p.draw()

    p = ReadAlign(fig, [w, w, w, w])
    p.delete(30, 70)
    p.draw()

    p = ReadAlign(fig, [w, 0, w, w])
    p.duplicate(30, 70)
    p.draw()

    p = OpticalMapAlign(fig, [2 * w, 0, w, w])
    p.draw()

    normalize_axes(root)

    image_name = mode + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
