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

from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.bed import Bed
from jcvi.formats.base import LineFile
from jcvi.apps.base import debug
from jcvi.utils.iter import pairwise

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
        self.label = args[2]
        self.bed = Bed(args[3])
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

    def __init__(self, ax, t, draw=True):

        # Copy the data from LayoutLine
        self.y = t.y
        self.sizes = sizes = t.sizes
        self.label = t.label
        self.color = t.color
        self.seqids = t.seqids
        self.bed = t.bed
        self.order = t.order
        self.order_in_chr = t.order_in_chr
        self.ax = ax

        self.xstart = xstart = .1
        self.gap = gap = .01
        rpad = .05
        span = 1 - xstart - rpad - gap * (len(sizes) - 1)
        total = sum(sizes.values())
        ratio = span / total

        self.ratio = ratio
        self.update_offsets()

        if draw:
            self.draw()

    def __str__(self):
        return self.label

    def draw(self):
        y = self.y
        color = self.color
        ax = self.ax
        xs = xstart = self.xstart
        gap = self.gap
        for sid in self.seqids:
            size = self.sizes[sid]
            xend = xstart + self.ratio * size
            hc = HorizontalChromosome(ax, xstart, xend, y, height=.02, fc=color)
            si = "".join(x for x in sid if x not in string.letters)
            si = str(int(si))
            xx = (xstart + xend) / 2
            ax.text(xx, y - .015, _(si), ha="center", va="top", color=color)
            xstart = xend + gap

        ax.text(xs / 2, y + gap, self.label, ha="center", color=color)

    def update_offsets(self):
        self.offsets = {}
        xs = self.xstart
        gap = self.gap
        for sid in self.seqids:
            size = self.sizes[sid]
            self.offsets[sid] = xs
            xs += self.ratio * size + gap

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


class PermutationSolver (object):
    """
    A heuristic method that rotates through the tracks, reorder the seqids
    according the seqids of another track - so as to have as few crossing shades
    as possible.
    """
    def __init__(self, tracks, layout):
        allblocks = {}
        for i, j, blocks in layout.edges:
            allblocks[(i, j)] = blocks
            a, b, c, d, score, orientation = zip(*blocks)
            blocks = zip(c, d, a, b, score, orientation)
            allblocks[(j, i)] = blocks

        logging.debug("Automatically resolve seqids ..")
        ntracks = len(tracks)
        # e.g. [2, 1, 0, 1, 2...], oscillating between 0 and ntracks
        pecking = range(ntracks - 1, 0, -1) + range(ntracks)[:-1]
        pecking *= 10  # How many cycles to go through
        score = 1e100  # Minimization problem
        for i, j in pairwise(pecking):
            blocks = allblocks[(i, j)]
            # Assume seqids in i already sorted, then sort seqids in j
            updated, score = self.sort_blocks(tracks[i], tracks[j],
                                              blocks, initialscore=score)
            if not updated:
                break

    def sort_blocks(self, atrack, btrack, blocks, initialscore=1e100):
        bseqids = set(btrack.seqids)
        finalscore = 0
        finalorder = []
        xstart = btrack.xstart
        gap = btrack.gap
        ratio = btrack.ratio
        order = btrack.order_in_chr
        sizes = btrack.sizes
        while bseqids:
            scores = defaultdict(int)
            for a, b, c, d, score, orientation in blocks:
                p = (atrack.get_coords(a)[0] + atrack.get_coords(b)[0]) / 2
                bseqid, c, f = order[c]
                bseqid, d, f = order[d]
                q = xstart + (c + d) / 2 * ratio
                scores[bseqid] += int(abs(p - q) * score) * 1000
            # Normalize the score over size
            for bseqid, score in scores.items():
                size = sizes[bseqid]
                scores[bseqid] /= size
            # Choose the best
            mv, mk = min((v, k) for k, v in scores.items() if k in bseqids)
            finalscore += mv
            finalorder.append(mk)
            bseqids.remove(mk)
            size = sizes[mk]
            xstart += size * ratio + gap

        logging.debug("Order {0} - {1}: initial score={2}, final score={3}".\
                      format(atrack, btrack, initialscore, finalscore))
        if finalscore < initialscore:
            btrack.seqids = finalorder
            btrack.update_offsets()
            return True, finalscore

        return False, -1


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
        tr = Track(root, lo, draw=False)
        tracks.append(tr)

    PermutationSolver(tracks, layout)
    ShadeManager(root, tracks, layout)

    for tr in tracks:
        tr.draw()  # this time for real

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
