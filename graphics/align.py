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

from jcvi.graphics.chromosome import Chromosome, HorizontalChromosome
from jcvi.graphics.base import plt, savefig, normalize_axes
from jcvi.apps.base import OptionParser


class BaseAlign (object):

    def __init__(self, fig, xywh, xmax=100):
        x, y, w, h = xywh
        self.ax = fig.add_axes([x, y, w, h])
        self.sax = fig.add_axes([x + .15 * w, y + .15 * h, w * .7, h * .7])
        self.amax = self.bmax = xmax
        self.a = [(0, xmax)]
        self.b = [(0, xmax)]
        self.apatch = self.bpatch = None

    def convert(self, pos, xmax):
        return .15 + pos * .7 / xmax

    def invert(self, a, b):
        self.a = [(0, a), (a, b), (b, self.amax)]
        self.b = [(0, a), (b, a), (b, self.bmax)]
        self.apatch = (self.convert(a, self.amax),
                       self.convert(b, self.amax))
        self.bpatch = (self.convert(a, self.bmax),
                       self.convert(b, self.bmax))

    def delete(self, a, b):
        self.bmax -= b - a
        self.a = [(0, a), (b, self.amax)]
        self.b = [(0, a), (a, self.bmax)]
        self.apatch = (self.convert(a, self.amax),
                       self.convert(b, self.amax))

    def duplicate(self, a, b, gap=0):
        self.bmax += b - a + gap
        self.a = [(0, b), (a, self.amax)]
        self.b = [(0, b), (b + gap, self.bmax)]
        self.apatch = (self.convert(a, self.amax),
                       self.convert(b, self.amax))
        self.bpatch = (self.convert(a, self.bmax),
                       self.convert(b, self.bmax),
                       self.convert(b + gap, self.bmax),
                       self.convert(2 * b - a + gap, self.bmax))


class PairwiseAlign (BaseAlign):

    def __init__(self, fig, xywh):
        super(PairwiseAlign, self).__init__(fig, xywh)

    def draw(self, width=.03):
        HorizontalChromosome(self.ax, .15, .85, .1, height=width,
                             patch=self.apatch, lw=2)
        Chromosome(self.ax, .1, .15, .85, width=width,
                   patch=self.bpatch, lw=2)
        for a, b in zip(self.a, self.b):
            self.sax.plot(a, b, "-", color="darkslategrey", lw=2)
        self.sax.set_xticklabels([])
        self.sax.set_yticklabels([])


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = p.set_image_options()

    if len(args) != 1:
        sys.exit(not p.print_help())

    mode, = args
    assert mode in ("dotplot", "reads", "om")

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    p = PairwiseAlign(fig, [0, 0, 1, 1])
    p.duplicate(20, 30, gap=10)
    p.draw()

    normalize_axes(root)

    image_name = mode + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
