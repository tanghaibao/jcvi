#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [dotplot|reads|om]

Illustrate three different types of alignments.
- Pairwise sequence alignment, aka, "dot plot"
- Read alignment, similar to the visualization of a BAM file
- Optical map alignment, matchings between restriction fragments
"""


import sys

from random import randint

from jcvi.utils.range import range_overlap
from jcvi.graphics.chromosome import Chromosome, HorizontalChromosome
from jcvi.graphics.glyph import GeneGlyph
from jcvi.graphics.base import plt, savefig, normalize_axes
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
                             self.ypad - .05, height=width,
                             patch=self.apatch, lw=2)
        Chromosome(self.ax, self.xpad - .05, self.ypad, 1 - self.ypad,
                   width=width, patch=self.bpatch, lw=2)
        for a, b in zip(self.a, self.b):
            self.sax.plot(a, b, "-", color="darkslategrey", lw=2)
        self.sax.set_xticklabels([])
        self.sax.set_yticklabels([])
        self.sax.set_xlim((1, self.amax))
        self.sax.set_ylim((1, self.bmax))


class ReadAlign (BaseAlign):

    def __init__(self, fig, xywh, xpad=.1, ypad=.25, readlen=6, gap=4):
        super(ReadAlign, self).__init__(fig, xywh, xpad, ypad)
        self.readlen = readlen
        self.gap = gap
        self.reads = []
        self.ymax = 12
        self.ntracks = 0
        self.layout(1, self.amax)

    def layout(self, start, end, maxtracks=6):
        readrange = 2 * self.readlen + self.gap
        end -= readrange
        assert start < end, "end must be > start + readlen"
        reads = []
        for x in xrange(100):
            pos = randint(start, end)
            reads.append(PairedRead(pos, readlen=self.readlen, gap=self.gap))

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
            r.set_y(self.ntracks + mi)
        reads = [x for x in reads if x.y is not None]
        self.reads += reads
        self.ntracks += len(track_ends)

    def remove(self, a, b, maxtracks=0):
        self.reads = [r for r in self.reads \
                      if not (range_overlap((0, a, b), (0, r.start, r.end)) \
                              and r.y >= maxtracks)]

    def draw(self, width=.03):
        HorizontalChromosome(self.ax, self.xpad, 1 - self.xpad,
                             self.ypad - width / 2, height=width,
                             patch=self.apatch, lw=2)
        for r in self.reads:
            r.draw(self.sax)
        normalize_axes(self.sax)
        self.sax.set_xlim((1, self.amax))
        self.sax.set_ylim((-1, self.ymax))

    def highlight(self, a, b):
        self.sax.plot((a, a), (-1, self.ntracks), "g-", lw=2)
        self.sax.plot((b, b), (-1, self.ntracks), "g-", lw=2)
        for r in self.reads:
            if range_overlap((0, a, b), (0, r.start, r.end)):
                r.set_color('r')

    def duplicate(self, a, b):
        self.apatch = (self.convert(a, self.amax),
                       self.convert(b, self.amax))
        self.layout(1, self.amax)
        self.remove(1, a, maxtracks=6)
        self.remove(b, self.amax, maxtracks=6)
        self.highlight(a, b)


class SingleRead (object):

    def __init__(self, start, readlen, sign=1):
        self.x1 = start
        self.x2 = start + sign * readlen
        self.y = None
        self.start, self.end = min(self.x1, self.x2), max(self.x1, self.x2)
        self.span = self.end - self.start + 1
        self.color = 'k'

    def draw(self, ax, height=.5):
        GeneGlyph(ax, self.x1, self.x2, self.y, height, tip=1,
                  color=self.color, gradient=True)

    def breakpoint(self, a, b):
        if self.start < a < self.end:
            self.end = a
            self.color = 'r'
        elif self.start < b < self.end:
            self.start = b
            self.color = 'r'


class PairedRead (object):

    def __init__(self, start, readlen, gap):
        self.r1 = SingleRead(start, readlen)
        i1 = start + readlen
        i2 = i1 + gap
        i3 = i2 + readlen
        self.r2 = SingleRead(i3, readlen, sign=-1)
        self.i1, self.i2 = i1, i2
        self.start = min(self.r1.start, self.r2.start)
        self.end = max(self.r1.end, self.r2.end)
        self.y = None

    def set_color(self, color):
        self.r1.color = self.r2.color = color

    def set_y(self, y):
        self.y = y
        self.r1.y = self.r2.y = y

    def draw(self, ax):
        self.r1.draw(ax)
        self.r2.draw(ax)
        ax.plot((self.i1, self.i2), (self.y, self.y), "-",
                 color="lightslategrey", lw=2)


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = p.set_image_options()

    if len(args) != 1:
        sys.exit(not p.print_help())

    mode, = args
    assert mode in ("dotplot", "reads", "om")

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    if mode == "dotplot":
        p = PairwiseAlign(fig, [0, 0, 1, 1])
        p.duplicate(30, 70, gap=5)
        p.draw()
    elif mode == "reads":
        p = ReadAlign(fig, [0, 0, 1, 1])
        p.duplicate(30, 70)
        p.draw()

    normalize_axes(root)

    image_name = mode + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
