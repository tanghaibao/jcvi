#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog seqids layout

Illustrate macrosynteny between tracks which represent individual genomes.

seqids contain the chromosomes to plot. Each line correspond to a track.
layout provides configuration for placement of tracks and mapping file between tracks.

Layout file example - first section specify how to draw each track. Then the "edges"
section specify which connections to draw.

# y, xstart, xend, rotation, color, label, va, bed, label_va
.6, .1, .4, 0, m, Grape, top, grape.bed, center
.4, .3, .6, 60, k, Athaliana, top, athaliana.bed, center
# edges
e, 0, 1, athaliana.grape.4x1.simple
"""


import sys
import logging

from jcvi.apps.base import OptionParser
from jcvi.compara.synteny import SimpleFile
from jcvi.formats.bed import Bed
from jcvi.graphics.chromosome import HorizontalChromosome
from jcvi.graphics.glyph import TextCircle
from jcvi.graphics.synteny import Shade
from jcvi.graphics.base import mpl, plt, savefig, markup, AbstractLayout


class LayoutLine(object):
    def __init__(self, row, delimiter=",", generank=True):
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]

        self.empty = False
        if len(args) < 8:
            self.empty = True
            return
        self.y = float(args[0])
        self.xstart = float(args[1])
        self.xend = float(args[2])
        self.rotation = int(args[3])
        self.color = args[4]
        self.label = args[5]
        self.va = args[6]
        self.bed = Bed(args[7])
        if len(args) == 9:
            self.label_va = args[8]
        else:
            self.label_va = "center"
        self.order = self.bed.order
        self.order_in_chr = self.bed.order_in_chr if generank else self.bed.bp_in_chr


class Layout(AbstractLayout):
    def __init__(self, filename, delimiter=",", generank=False):
        super(Layout, self).__init__(filename)
        fp = open(filename)
        self.edges = []
        for row in fp:
            if row[0] == "#":
                continue
            if row[0] == "e":
                args = row.rstrip().split(delimiter)
                args = [x.strip() for x in args]
                i, j, fn = args[1:4]
                if len(args) == 5 and args[4]:
                    samearc = args[4]
                else:
                    samearc = "below"
                i, j = int(i), int(j)
                assert args[0] == "e"
                blocks = self.parse_blocks(fn, i)
                self.edges.append((i, j, blocks, samearc))
            else:
                self.append(LayoutLine(row, delimiter=delimiter, generank=generank))

        self.assign_colors()

    def parse_blocks(self, simplefile, i):
        order = self[i].order
        return SimpleFile(simplefile, order=order).blocks


MaxSeqids = 16  # above which no labels are written


def make_circle_name(sid, rev):
    """Extract a succinct label based on sid.

    If there are numbers to be extracted, returns the first group of number.
    Otherwise, the first letter is returned.

    If sid is in rev, then '-' gets appended to the label.

    Args:
        sid (str): seqid
        rev (set[str]): Set of seqids that are reversed

    Returns:
        str: Single letter label for the sid
    """
    import re

    in_reverse = sid in rev
    sid = sid.rsplit("_", 1)[-1]
    si = re.findall(r"\d+", sid)
    if si:
        si = str(int(si[0]))
    else:
        si = sid[0].upper()
    if in_reverse:
        si += "-"
    return si


class Track(object):
    def __init__(
        self,
        ax,
        t,
        gap=0.01,
        height=0.01,
        lw=1,
        draw=True,
        roundrect=False,
    ):

        self.empty = t.empty
        if t.empty:
            return

        # Copy the data from LayoutLine
        self.y = t.y
        self.sizes = sizes = t.sizes
        self.label = t.label
        self.rotation = t.rotation
        self.va = t.va
        self.label_va = t.label_va
        self.color = t.color if t.color != "None" else None
        self.seqids = t.seqids
        self.bed = t.bed
        self.order = t.order
        self.order_in_chr = t.order_in_chr
        self.rev = t.rev
        self.ax = ax
        self.height = height

        self.xstart = xstart = t.xstart
        self.xend = t.xend

        # Rotation transform
        self.x = x = (self.xstart + self.xend) / 2
        y = self.y
        self.tr = (
            mpl.transforms.Affine2D().rotate_deg_around(x, y, self.rotation)
            + ax.transAxes
        )
        self.inv = ax.transAxes.inverted()

        nseqids = len(self.seqids)
        if nseqids > MaxSeqids:
            gap = min(gap, gap * MaxSeqids / nseqids + 0.001)
        self.gap = gap

        rpad = 1 - t.xend
        span = 1 - xstart - rpad - gap * (len(sizes) - 1)
        self.total = total = sum(sizes.values())
        ratio = span / total

        self.ratio = ratio
        self.update_offsets()
        self.lw = lw

        if draw:
            self.draw(roundrect=roundrect)

    def __str__(self):
        return self.label

    def draw(
        self, roundrect=False, plot_label=True, plot_circles=True, pad=0.03, vpad=0.09
    ):
        if self.empty:
            return

        y = self.y
        color = self.color
        ax = self.ax
        xstart = self.xstart
        gap = self.gap
        va = self.va
        nseqids = len(self.seqids)
        tr = self.tr

        for i, sid in enumerate(self.seqids):
            size = self.sizes[sid]
            rsize = self.ratio * size
            xend = xstart + rsize
            hc = HorizontalChromosome(
                ax,
                xstart,
                xend,
                y,
                height=self.height,
                lw=self.lw,
                fc=color,
                roundrect=roundrect,
            )
            hc.set_transform(tr)
            si = make_circle_name(sid, self.rev)
            xx = (xstart + xend) / 2
            xstart = xend + gap

            step = 2 if nseqids <= 40 else 10
            if nseqids >= 2 * MaxSeqids and (i + 1) % step != 0:
                continue
            if nseqids < 5:
                continue

            hpad = -pad if va == "bottom" else pad
            if plot_circles:
                TextCircle(
                    ax,
                    xx,
                    y + hpad,
                    si,
                    fc="w",
                    color=color,
                    size=10,
                    transform=tr,
                )

        label = markup(self.label)
        c = color if color != "gainsboro" else "k"
        if plot_label:
            if self.label_va == "top":
                x, y = self.x, self.y + vpad
            elif self.label_va == "bottom":
                x, y = self.x, self.y - vpad
            else:  # "center"
                x, y = self.xstart - vpad / 2, self.y
            ax.text(x, y, label, ha="center", va="center", color=c, transform=tr)

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
        seqid, i, _ = order_in_chr[gene]
        if seqid not in self.offsets:
            return [None, None]

        x = self.offsets[seqid]
        if seqid in self.rev:
            x += self.ratio * (self.sizes[seqid] - i - 1)
        else:
            x += self.ratio * i
        y = self.y
        x, y = self.tr.transform((x, y))
        x, y = self.inv.transform((x, y))

        return [x, y]


class ShadeManager(object):
    def __init__(self, ax, tracks, layout, heightpad=0, style="curve"):
        self.style = style
        for i, j, blocks, samearc in layout.edges:
            # if same track (duplication shades), shall we draw above or below?
            # samearc = "above" if i == j and i == 0 else "below"
            self.draw_blocks(
                ax, blocks, tracks[i], tracks[j], samearc=samearc, heightpad=heightpad
            )

    def draw_blocks(self, ax, blocks, atrack, btrack, samearc="below", heightpad=0):
        for a, b, c, d, _, _, highlight in blocks:
            p = atrack.get_coords(a), atrack.get_coords(b)
            q = btrack.get_coords(c), btrack.get_coords(d)
            if p[0] is None or q[0] is None:
                continue

            ymid = (atrack.y + btrack.y) / 2
            px, qx = p[0][0], q[0][0]
            xdist = abs(px - qx) if px and qx else 0.5
            pad = 0.09 * xdist / 0.5
            if atrack.y == btrack.y:
                if samearc == "below":
                    ymid = atrack.y - pad
                else:
                    ymid = atrack.y + pad
            if heightpad:
                if atrack.y < btrack.y:
                    p[0][1] = p[1][1] = atrack.y + heightpad
                    q[0][1] = q[1][1] = btrack.y - heightpad
                else:
                    p[0][1] = p[1][1] = atrack.y - heightpad
                    q[0][1] = q[1][1] = btrack.y + heightpad

            zorder = 2 if highlight else 1
            lw = 1 if highlight else 0
            Shade(
                ax,
                p,
                q,
                ymid,
                highlight=highlight,
                alpha=1,
                fc="gainsboro",
                ec="gainsboro",
                lw=lw,
                zorder=zorder,
                style=self.style,
            )


class Karyotype(object):
    def __init__(
        self,
        fig,
        root,
        seqidsfile,
        layoutfile,
        gap=0.01,
        height=0.01,
        lw=1,
        generank=True,
        sizes=None,
        heightpad=0,
        roundrect=False,
        plot_label=True,
        plot_circles=True,
        shadestyle="curve",
    ):
        layout = Layout(layoutfile, generank=generank)

        fp = open(seqidsfile)
        # Strip the reverse orientation tag for e.g. chr3-
        di = lambda x: x[:-1] if x[-1] == "-" else x
        for i, row in enumerate(fp):
            if row[0] == "#":
                continue
            t = layout[i]
            # There can be comments in seqids file:
            # https://github.com/tanghaibao/jcvi/issues/335
            seqids = row.split("#", 1)[0].rstrip().split(",")
            t.rev = set(x[:-1] for x in seqids if x[-1] == "-")
            seqids = [di(x) for x in seqids]
            if t.empty:
                continue

            bed = t.bed
            self.generank = generank
            if generank:
                sz = dict((x, len(list(bed.sub_bed(x)))) for x in seqids)
            else:
                sz = sizes or dict(
                    (x, max(z.end for z in list(bed.sub_bed(x)))) for x in seqids
                )
                assert sz is not None, "sizes not available and cannot be inferred"
            t.seqids = seqids
            # validate if all seqids are non-empty
            for k, v in sz.items():
                if v == 0:
                    logging.error("Size of `%s` is empty. Please check", k)
            t.sizes = sz

        tracks = []
        for lo in layout:
            if lo.empty:
                continue
            tr = Track(root, lo, gap=gap, height=height, lw=lw, draw=False)
            tracks.append(tr)

        ShadeManager(root, tracks, layout, heightpad=heightpad, style=shadestyle)

        for tr in tracks:
            tr.draw(
                roundrect=roundrect, plot_label=plot_label, plot_circles=plot_circles
            )

        self.tracks = tracks
        self.layout = layout


def main():
    p = OptionParser(__doc__)
    p.add_option(
        "--basepair",
        default=False,
        action="store_true",
        help="Use base pair position instead of gene rank",
    )
    p.add_option(
        "--nocircles",
        default=False,
        action="store_true",
        help="Do not plot chromosome circles",
    )
    p.add_option(
        "--shadestyle",
        default="curve",
        choices=Shade.Styles,
        help="Style of syntenic wedges",
    )
    opts, args, iopts = p.set_image_options(figsize="8x7")

    if len(args) != 2:
        sys.exit(not p.print_help())

    seqidsfile, layoutfile = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Karyotype(
        fig,
        root,
        seqidsfile,
        layoutfile,
        plot_circles=(not opts.nocircles),
        shadestyle=opts.shadestyle,
        generank=(not opts.basepair),
    )

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "karyotype"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main()
