#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog anchorsfile all.bed seqids layout

Illustrate macrosynteny between tracks which represent individual genomes.

seqids contain the chromosomes to plot. Each line correspond to a track.
layout provides configuration for placement of tracks and mapping file between tracks.
"""


import sys
import string
import os.path as op
import logging

from optparse import OptionParser

from jcvi.graphics.chromosome import HorizontalChromosome
from jcvi.graphics.base import plt, _, set_image_options
from jcvi.formats.bed import Bed
from jcvi.formats.base import LineFile
from jcvi.apps.base import debug
debug()


class LayoutLine (object):

    def __init__(self, row, delimiter=','):
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.y = float(args[0])
        self.color = args[1]


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
                a, b, fn = args[1:4]
                a, b = int(a), int(b)
                assert args[0] == 'e'
                blocks = self.parse_blocks(fn)
                self.edges.append((a, b, blocks))
            else:
                self.append(LayoutLine(row, delimiter=delimiter))

    def parse_blocks(self, simplefile):
        fp = open(simplefile)
        blocks = []
        for row in fp:
            a, b, c, d, score, orientation = row.split()
            score = int(score)
            blocks.append((a, b, c, d, score, orientation))


class Track (object):

    def __init__(self, ax, y, seqids, sizes, color="w"):

        assert len(seqids) == len(sizes)

        xstart = .05
        gap = .04
        span = .9 - gap * (len(sizes) - 1)
        total = sum(sizes.values())
        ratio = span / total

        for sid in seqids:
            size = sizes[sid]
            xend = xstart + ratio * size
            hc = HorizontalChromosome(ax, xstart, xend, y, height=.02, fc=color)
            si = "".join(x for x in sid if x not in string.letters)
            si = str(int(si))
            xx = (xstart + xend) / 2
            ax.text(xx, y - .015, _(si), ha="center", va="top", color=color)
            xstart = xend + gap

def main():
    p = OptionParser(__doc__)
    opts, args, iopts = set_image_options(p, figsize="8x4")

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, seqidsfile, layoutfile = args
    bed = Bed(bedfile)
    lo = Layout(layoutfile)
    pf = bedfile.rsplit(".", 1)[0]

    fp = open(seqidsfile)
    tseqids = []
    tsizes = []
    for row in fp:
        seqids = row.rstrip().split(",")
        sizes = dict((x, len(list(bed.sub_bed(x)))) for x in seqids)
        tseqids.append(seqids)
        tsizes.append(sizes)

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    ntracks = len(lo)
    for t, seqids, sizes in zip(lo, tseqids, tsizes):
        Track(root, t.y, seqids, sizes, color=t.color)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


if __name__ == '__main__':
    main()
