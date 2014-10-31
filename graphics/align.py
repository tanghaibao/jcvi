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

from jcvi.graphics.base import plt, savefig, normalize_axes
from jcvi.apps.base import OptionParser


class BaseAlign (object):

    def __init__(self, ax, xmax=100):
        self.ax = ax
        self.amax = self.bmax = xmax
        self.a = [(0, xmax)]
        self.b = [(0, xmax)]

    def invert(self, a, b):
        self.a = [(0, a), (a, b), (b, self.amax)]
        self.b = [(0, a), (b, a), (b, self.bmax)]

    def delete(self, a, b):
        self.bmax -= b - a
        self.a = [(0, a), (b, self.amax)]
        self.b = [(0, a), (a, self.bmax)]

    def duplicate(self, a, b):
        self.bmax += b - a
        self.a = [(0, b), (a, self.amax)]
        self.b = [(0, b), (b, self.bmax)]


class PairwiseAlign (BaseAlign):

    def __init__(self, ax):
        super(PairwiseAlign, self).__init__(ax)

    def draw(self):
        ax = self.ax
        for a, b in zip(self.a, self.b):
            ax.plot(a, b, "-", color="darkslategrey")


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = p.set_image_options()

    if len(args) != 1:
        sys.exit(not p.print_help())

    mode, = args
    assert mode in ("dotplot", "reads", "om")

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    ax = fig.add_axes([.1, .1, .8, .8])
    p = PairwiseAlign(ax)
    p.invert(20, 30)
    p.draw()

    normalize_axes(root)

    image_name = mode + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
