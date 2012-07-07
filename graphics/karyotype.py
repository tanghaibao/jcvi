#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog seqids layout

Illustrate macrosynteny between tracks which represent individual genomes.

seqids contain the chromosomes to plot. Each line correspond to a track.
layout provides configuration for placement of tracks and mapping file between tracks.
"""


import sys
import string
import os.path as op
import logging

from optparse import OptionParser

from jcvi.formats.bed import Bed
from jcvi.formats.base import LineFile
from jcvi.apps.base import debug

from jcvi.graphics.chromosome import HorizontalChromosome
from jcvi.graphics.synteny import Shade
from jcvi.graphics.base import plt, _, set_image_options
debug()


class LayoutLine (object):

    def __init__(self, row, delimiter=','):
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.y = float(args[0])
        self.color = args[1]
        self.bed = Bed(args[2])
        self.order = self.bed.order
        self.order_in_chr = self.bed.order_in_chr


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
                i, j, fn = args[1:4]
                i, j = int(i), int(j)
                assert args[0] == 'e'
                blocks = self.parse_blocks(fn, i)
                self.edges.append((i, j, blocks))
            else:
                self.append(LayoutLine(row, delimiter=delimiter))

    def parse_blocks(self, simplefile, i):
        order = self[i].order
        # Sometimes the simplefile has query and subject wrong
        fp = open(simplefile)
        blocks = []
        for row in fp:
            a, b, c, d, score, orientation = row.split()
            if a not in order:
                a, b, c, d = c, d, a, b
            if orientation == '-':
                c, d = d, c
            score = int(score)
            blocks.append((a, b, c, d, score, orientation))
        return blocks


class Track (object):

    def __init__(self, ax, t):

        # Copy the data from LayoutLine
        self.y = y = t.y
        seqids = t.seqids
        sizes = t.sizes
        self.bed = t.bed
        self.order = t.order
        self.order_in_chr = t.order_in_chr
        color = t.color

        xstart = .05
        gap = .01
        span = .9 - gap * (len(sizes) - 1)
        total = sum(sizes.values())
        ratio = span / total
        offsets = {}

        for sid in seqids:
            size = sizes[sid]
            offsets[sid] = xstart
            xend = xstart + ratio * size
            hc = HorizontalChromosome(ax, xstart, xend, y, \
                                      height=.02, fc=color)
            si = "".join(x for x in sid if x not in string.letters)
            si = str(int(si))
            xx = (xstart + xend) / 2
            ax.text(xx, y - .015, _(si), ha="center", va="top", color=color)
            xstart = xend + gap

        self.offsets = offsets
        self.ratio = ratio

    def get_coords(self, gene):
        order_in_chr = self.order_in_chr
        seqid, i, f = order_in_chr[gene]
        x = self.offsets[seqid] + self.ratio * i
        y = self.y
        return x, y


class ShadeManager (object):

    def __init__(self, ax, tracks, layout):
        for i, j, blocks in layout.edges:
            self.draw_blocks(ax, blocks, tracks[i], tracks[j])

    def draw_blocks(self, ax, blocks, atrack, btrack):
        for a, b, c, d, score, orientation in blocks:
            p = atrack.get_coords(a), atrack.get_coords(b)
            q = btrack.get_coords(c), btrack.get_coords(d)
            ymid = (atrack.y + btrack.y) / 2

            Shade(ax, p, q, ymid, highlight=False)


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = set_image_options(p, figsize="8x4")

    if len(args) != 2:
        sys.exit(not p.print_help())

    seqidsfile, layoutfile = args
    layout = Layout(layoutfile)

    fp = open(seqidsfile)
    for i, row in enumerate(fp):
        t = layout[i]
        seqids = row.rstrip().split(",")
        bed = t.bed
        sizes = dict((x, len(list(bed.sub_bed(x)))) for x in seqids)
        t.seqids = seqids
        t.sizes = sizes

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    ntracks = len(layout)
    tracks = []
    for lo in layout:
        tr = Track(root, lo)
        tracks.append(tr)

    ShadeManager(root, tracks, layout)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "out"
    image_name = pf + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


if __name__ == '__main__':
    main()
