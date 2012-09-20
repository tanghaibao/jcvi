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
import random

from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.bed import Bed
from jcvi.formats.base import LineFile
from jcvi.apps.base import debug
from jcvi.utils.iter import pairwise

from jcvi.graphics.chromosome import HorizontalChromosome
from jcvi.graphics.glyph import TextCircle
from jcvi.graphics.synteny import Shade
from jcvi.graphics.base import plt, _, set_image_options, savefig
debug()


class LayoutLine (object):

    def __init__(self, row, delimiter=','):
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.y = float(args[0])
        self.xstart = float(args[1])
        self.xend = float(args[2])
        self.color = args[3]
        self.label = args[4]
        self.va = args[5]
        self.bed = Bed(args[6])
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
            hl = ("*" in row)
            if hl:
                hl, row = row.split("*", 1)
                hl = hl or "r"
            a, b, c, d, score, orientation = row.split()
            if a not in order:
                a, b, c, d = c, d, a, b
            if orientation == '-':
                c, d = d, c
            score = int(score)
            blocks.append((a, b, c, d, score, orientation, hl))
        return blocks


class Track (object):

    def __init__(self, ax, t, draw=True):

        # Copy the data from LayoutLine
        self.y = t.y
        self.sizes = sizes = t.sizes
        self.label = t.label
        self.va = t.va
        self.color = t.color
        self.seqids = t.seqids
        self.bed = t.bed
        self.order = t.order
        self.order_in_chr = t.order_in_chr
        self.ax = ax

        self.xstart = xstart = t.xstart
        self.xend = t.xend
        gap = .01
        nseqids = len(self.seqids)
        if nseqids > 40:
            gap /= (nseqids / 40.)
        self.gap = gap

        rpad = 1 - t.xend
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
        va = self.va
        nseqids = len(self.seqids)
        for i, sid in enumerate(self.seqids):
            size = self.sizes[sid]
            rsize = self.ratio * size
            xend = xstart + rsize
            hc = HorizontalChromosome(ax, xstart, xend, y, height=.01, fc=color)
            sid = sid.rsplit("_", 1)[-1]
            si = "".join(x for x in sid if x not in string.letters)
            si = str(int(si))
            xx = (xstart + xend) / 2
            xstart = xend + gap

            if nseqids > 40:
                continue

            pad = .02
            if va == "bottom":
                pad = - pad
            TextCircle(ax, xx, y + pad, _(si), radius=.01,
                       fc="w", color=color, size=10)

        xp = self.xstart / 2 if self.xstart <= .5 else (1 + self.xend) / 2
        ax.text(xp, y + gap, self.label, ha="center", color=color)

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
        for a, b, c, d, score, orientation, highlight in blocks:
            p = atrack.get_coords(a), atrack.get_coords(b)
            q = btrack.get_coords(c), btrack.get_coords(d)
            ymid = (atrack.y + btrack.y) / 2
            if atrack.y == btrack.y:
                ymid = atrack.y - .09

            zorder = 2 if highlight else 1
            Shade(ax, p, q, ymid, highlight=highlight, alpha=1, fc="gainsboro",
                        lw=0, zorder=zorder)


class PermutationSolver (object):
    """
    A heuristic method that rotates through the tracks, reorder the seqids
    according the seqids of another track - so as to have as few crossing shades
    as possible.
    """
    def __init__(self, tracks, layout):
        self.tracks = tracks
        allblocks = {}
        for i, j, blocks in layout.edges:
            allblocks[(i, j)] = blocks
            a, b, c, d, score, orientation, highlight = zip(*blocks)
            blocks = zip(c, d, a, b, score, orientation, highlight)
            allblocks[(j, i)] = blocks
        self.allblocks = allblocks

        logging.debug("Automatically resolve seqids ..")
        ntracks = len(tracks)
        # e.g. [2, 1, 0, 1, 2...], oscillating between 0 and ntracks
        score = 1e100  # Minimization problem
        ff = range(ntracks)
        rr = ff[::-1]
        peckings = [ff, rr] * 20  # Max 2 x 100 iterations in total
        for pecking in peckings:
            updated, score = self.forward(pecking, initialscore=score)
            if not updated:
                break

    def forward(self, pecking, initialscore=1e100):
        # Save current seqids in case this run fails
        backups = [x.seqids for x in self.tracks]
        finalscore = 0
        allblocks = self.allblocks
        tracks = self.tracks
        opecking = list(pairwise(pecking))
        #random.shuffle(opecking)
        for i, j in opecking:
            blocks = allblocks[(i, j)]
            # Assume seqids in i already sorted, then sort seqids in j
            score = self.sort_blocks(tracks[i], tracks[j], blocks)
            finalscore += score

        updated = finalscore < initialscore
        label = "updated" if updated else "rejected"
        logging.debug("Initial score={0}, final score={1}, {2}".\
                      format(initialscore, finalscore, label))
        if finalscore < initialscore:
            return True, finalscore

        for backup, tr in zip(backups, self.tracks):
            tr.seqids = backup
            tr.update_offsets()
        return False, initialscore

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
            for a, b, c, d, score, orientation, highlight in blocks:
                p = (atrack.get_coords(a)[0] + atrack.get_coords(b)[0]) / 2
                bseqid, c, f = order[c]
                bseqid, d, f = order[d]
                q = xstart + (c + d) / 2 * ratio
                bscore = int(abs(p - q) * score * 1000)
                scores[bseqid] += bscore
            # Normalize the score over size
            for bseqid, score in scores.items():
                size = sizes[bseqid]
                scores[bseqid] /= size
            # Choose the best
            ss = [(v, k) for k, v in scores.items() if k in bseqids]
            if not ss:
                break
            mv, mk = min(ss)
            finalscore += mv
            finalorder.append(mk)
            bseqids.remove(mk)
            size = sizes[mk]
            xstart += size * ratio + gap

        btrack.seqids = finalorder
        btrack.update_offsets()

        return finalscore


class Karyotype (object):

    def __init__(self, fig, root, seqidsfile, layoutfile, auto=False):

        layout = Layout(layoutfile)

        fp = open(seqidsfile)
        for i, row in enumerate(fp):
            if row[0] == '#':
                continue
            t = layout[i]
            seqids = row.rstrip().split(",")
            bed = t.bed
            sizes = dict((x, len(list(bed.sub_bed(x)))) for x in seqids)
            t.seqids = seqids
            t.sizes = sizes

        tracks = []
        for lo in layout:
            tr = Track(root, lo, draw=False)
            tracks.append(tr)

        if auto:
            PermutationSolver(tracks, layout)
            autofile = seqidsfile + ".auto"
            fw = open(autofile, "w")
            for tr in tracks:
                print >> fw, ",".join(tr.seqids)
            fw.close()
            logging.debug("Auto seqids written to `{0}`.".format(autofile))

        ShadeManager(root, tracks, layout)

        for tr in tracks:
            tr.draw()  # this time for real

        self.tracks = tracks


def main():
    p = OptionParser(__doc__)
    p.add_option("--auto", action="store_true",
                 help="Automatically adjust seqids [default: %default]")
    opts, args, iopts = set_image_options(p, figsize="8x7")

    if len(args) != 2:
        sys.exit(not p.print_help())

    seqidsfile, layoutfile = args
    auto = opts.auto

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Karyotype(fig, root, seqidsfile, layoutfile, auto=auto)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "out"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
