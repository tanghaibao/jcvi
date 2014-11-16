#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Legacy script to plot distribution of certain classes onto chromosomes. Adapted
from the script used in the Tang et al. PNAS 2010 paper, sigma figure.
"""

import sys
import logging

from math import ceil
from itertools import groupby

import numpy as np

from jcvi.formats.base import DictFile
from jcvi.formats.bed import Bed
from jcvi.graphics.glyph import RoundRect
from jcvi.graphics.base import plt, Rectangle, Polygon, CirclePolygon, savefig
from jcvi.graphics.glyph import BaseGlyph, plot_cap
from jcvi.apps.base import OptionParser


class Chromosome (BaseGlyph):

    def __init__(self, ax, x, y1, y2, width=.015, ec="k", patch=None,
                 patchcolor='lightgrey', lw=1, fc="k", zorder=2):
        """
        Chromosome with positions given in (x, y1) => (x, y2)

        The chromosome can also be patched, e.g. to show scaffold composition in
        alternating shades. Use a list of starting locations to segment.
        """
        y1, y2 = sorted((y1, y2))
        super(Chromosome, self).__init__(ax)
        pts, r = self.get_pts(x, y1, y2, width)
        self.append(Polygon(pts, fill=False, lw=lw, ec=ec, zorder=zorder))
        if patch:
            rr = r * .9  # Shrink a bit for the patches
            for i in xrange(0, len(patch), 2):
                if i + 1 > len(patch) - 1:
                    continue
                p1, p2 = patch[i], patch[i + 1]
                self.append(Rectangle((x - rr, p1), 2 * rr, p2 - p1, lw=0,
                             fc=patchcolor))

        self.add_patches()

    def get_pts(self, x, y1, y2, width):
        w = width / 2
        r = width / (3 ** .5)

        pts = []
        pts += plot_cap((x, y1 + r), np.radians(range(210, 330)), r)
        pts += [[x + w, y1 + r / 2], [x + w, y2 - r / 2]]
        pts += plot_cap((x, y2 - r), np.radians(range(30, 150)), r)
        pts += [[x - w, y2 - r / 2], [x - w, y1 + r / 2]]

        return pts, r


class HorizontalChromosome (BaseGlyph):

    def __init__(self, ax, x1, x2, y, height=.015, ec="k", patch=None,
                 patchcolor='lightgrey', lw=1, fc=None,
                 zorder=2, roundrect=False):
        """
        Horizontal version of the Chromosome glyph above.
        """
        x1, x2 = sorted((x1, x2))
        super(HorizontalChromosome, self).__init__(ax)
        pts, r = self.get_pts(x1, x2, y, height)
        if roundrect:
            RoundRect(ax, (x1, y - height * .5), x2 - x1, height, fill=False,
                      lw=lw, ec=ec, zorder=zorder + 1)
        else:
            self.append(Polygon(pts, fill=False, lw=lw, ec=ec, zorder=zorder))

        if fc:
            pts, r = self.get_pts(x1, x2, y, height / 2)
            if roundrect:
                RoundRect(ax, (x1, y - height / 4), x2 - x1, height / 2, fc=fc,
                          lw=0, zorder=zorder)
            else:
                self.append(Polygon(pts, fc=fc, lw=0, zorder=zorder))
        if patch:
            rr = r * .9  # Shrink a bit for the patches
            for i in xrange(0, len(patch), 2):
                if i + 1 > len(patch) - 1:
                    continue
                p1, p2 = patch[i], patch[i + 1]
                self.append(Rectangle((p1, y - rr), p2 - p1, 2 * rr, lw=0,
                             fc=patchcolor))

        self.add_patches()

    def get_pts(self, x1, x2, y, height):
        h = height / 2
        r = height / (3 ** .5)

        if x2 - x1 < 2 * height:  # rectangle for small chromosomes
            return [[x1, y + h], [x1, y - h], [x2, y - h], [x2, y + h]], r

        pts = []
        pts += plot_cap((x1 + r, y), np.radians(range(120, 240)), r)
        pts += [[x1 + r / 2, y - h], [x2 - r / 2, y - h]]
        pts += plot_cap((x2 - r, y), np.radians(range(-60, 60)), r)
        pts += [[x2 - r / 2, y + h], [x1 + r / 2, y + h]]

        return pts, r


class ChromosomeWithCentromere (object):

    def __init__(self, ax, x, y1, y2, y3, width=.015, fc="k", fill=False, zorder=2):
        """
        Chromosome with centromeres at y2 position
        """
        pts = []
        r = width * .5
        pts += plot_cap((x, y1 - r), np.radians(range(180)), r)
        pts += [[x - r, y1 - r], [x - r, y2 + r]]
        pts += plot_cap((x, y2 + r), np.radians(range(180, 360)), r)
        pts += [[x + r, y2 + r], [x + r, y1 - r]]
        ax.add_patch(Polygon(pts, fc=fc, fill=fill, zorder=zorder))
        pts = []
        pts += plot_cap((x, y2 - r), np.radians(range(180)), r)
        pts += [[x - r, y2 - r], [x - r, y3 + r]]
        pts += plot_cap((x, y3 + r), np.radians(range(180, 360)), r)
        pts += [[x + r, y3 + r], [x + r, y2 - r]]
        ax.add_patch(Polygon(pts, fc=fc, fill=fill, zorder=zorder))
        ax.add_patch(CirclePolygon((x, y2), radius=r * .5,
            fc="k", ec="k", zorder=zorder))


class ChromosomeMap (object):
    """
    Line plots along the chromosome.
    """
    def __init__(self, fig, root, xstart, xend, ystart, yend, pad, ymin, ymax, bins,
                    title, subtitle, patchstart=None):

        width, height = xend - xstart, yend - ystart

        y = ystart - pad
        HorizontalChromosome(root, xstart, xend, y, patch=patchstart, height=.03)

        # Gauge
        lsg = "lightslategrey"
        root.plot([xstart - pad, xstart - pad], [ystart, ystart + height],
                    lw=2, color=lsg)
        root.plot([xend + pad, xend + pad], [ystart, ystart + height],
                    lw=2, color=lsg)
        root.text((xstart + xend) / 2, ystart + height + 2 * pad, title,
                    ha="center", va="center", color=lsg)

        iv = (ymax - ymin) / bins
        iv_height = height / bins
        val = ymin
        yy = ystart
        while val <= ymax:
            root.text(xstart - 2 * pad, yy, str(val), ha="right", va="center", size=10)
            val += iv
            yy += iv_height

        root.text((xstart + xend) / 2, y - .05, subtitle, ha="center", va="center", color=lsg)

        self.axes = fig.add_axes([xstart, ystart, width, height])


class GeneticMap (BaseGlyph):

    def __init__(self, ax, x, y1, y2, markers, unit="cM",
                 tip=.008, fc="k", flip=False):
        # tip = length of the ticks
        y1, y2 = sorted((y1, y2))
        ax.plot([x, x], [y1, y2], '-', color=fc, lw=2)
        max_marker_name, max_chr_len = max(markers, key=lambda x: x[-1])
        r = y2 - y1
        ratio = r / max_chr_len
        marker_pos = {}
        for marker_name, cm in markers:
            yy = (y1 + ratio * cm) if flip else (y2 - ratio * cm)
            ax.plot((x - tip, x + tip), (yy, yy), "-", color=fc)
            marker_pos[marker_name] = yy
        self.marker_pos = marker_pos

        t = tip / 2
        end_cm_labels = ((y2 + t, max_chr_len, "bottom"), (y1 - t, 0, "top")) if flip \
                        else ((y2 + t, 0, "bottom"), (y1 - t, max_chr_len, "top"))
        for yy, cm, va in end_cm_labels:
            label = "{0} {1}".format(int(cm), unit)
            ax.text(x, yy, label, color="gray", va=va, ha="center")


class Gauge (BaseGlyph):

    def __init__(self, ax, x, y1, y2, max_chr_len, step=1e6,
                 tip=.008, extra=.006, fc="b"):
        # tip = length of the ticks
        # extra = offset for the unit label
        ax.plot([x, x], [y1, y2], '-', color=fc, lw=2)
        r = y2 - y1
        yy = y2
        gauge = int(ceil(max_chr_len / step))
        ratio = r / max_chr_len
        yinterval = 2 * ratio * step
        for g in xrange(0, gauge, 2):
            if g % 10:
                ax.plot((x, x + tip), (yy, yy), "-", color=fc)
            else:
                ax.plot((x - tip, x + tip), (yy, yy), '-', color=fc, lw=2)
                ax.text(x + tip + extra, yy, g, color="gray", va="center")
            yy -= yinterval
        ax.text(x, yy - .03, "Mb", color="gray", va="center")


def canvas2px(coord, dmn, dpi):
    """
    Convert matplotlib canvas coordinate to pixels
    """
    return int(round(coord * dmn * dpi))


def write_ImageMapLine(tlx, tly, brx, bry, w, h, dpi, chr, segment_start, segment_end):
    """
    Write out an image map area line with the coordinates passed to this
    function
    <area shape="rect" coords="tlx,tly,brx,bry" href="#chr7" title="chr7:100001..500001">
    """
    tlx, brx = [canvas2px(x, w, dpi) for x in (tlx, brx)]
    tly, bry = [canvas2px(y, h, dpi) for y in (tly, bry)]
    chr, bac_list = chr.split(':')
    return '<area shape="rect" coords="' + \
           ",".join(str(x) for x in (tlx, tly, brx, bry)) \
           + '" href="#' + chr + '"' \
           + ' title="' + chr + ':' + str(segment_start) + '..' + str(segment_end) + '"' \
           + ' />'


def main():
    """
    %prog bedfile id_mappings

    Takes a bedfile that contains the coordinates of features to plot on the
    chromosomes, and `id_mappings` file that map the ids to certain class. Each
    class will get assigned a unique color. `id_mappings` file is optional (if
    omitted, will not paint the chromosome features, except the centromere).
    """
    p = OptionParser(main.__doc__)
    p.add_option("--title", default="Medicago truncatula v3.5",
            help="title of the image [default: `%default`]")
    p.add_option("--gauge", default=False, action="store_true",
            help="draw a gauge with size label [default: %default]")
    p.add_option("--imagemap", default=False, action="store_true",
            help="generate an HTML image map associated with the image [default: %default]")
    p.add_option("--winsize", default=50000, type="int",
            help="if drawing an imagemap, specify the window size (bases) of each map element "
                 "[default: %default bp]")
    p.add_option("--empty",
            help="Write legend for unpainted region")
    opts, args, iopts = p.set_image_options(figsize="6x6", dpi=300)

    if len(args) not in (1, 2):
        sys.exit(p.print_help())

    bedfile = args[0]
    mappingfile = None
    if len(args) == 2:
        mappingfile = args[1]

    winsize = opts.winsize
    imagemap = opts.imagemap
    w, h = iopts.w, iopts.h
    dpi = iopts.dpi

    prefix = bedfile.rsplit(".", 1)[0]
    figname = prefix + "." + opts.format
    if imagemap:
        imgmapfile = prefix + '.map'
        mapfh = open(imgmapfile, "w")
        print >> mapfh, '<map id="' + prefix + '">'

    if mappingfile:
        mappings = DictFile(mappingfile, delimiter="\t")
        classes = sorted(set(mappings.values()))
        logging.debug("A total of {0} classes found: {1}".format(len(classes),
            ','.join(classes)))
    else:
        mappings = {}
        classes = []
        logging.debug("No classes registered (no id_mappings given).")

    mycolors = "rgbymc"
    class_colors = dict(zip(classes, mycolors))

    bed = Bed(bedfile)
    chr_lens = {}
    centromeres = {}
    for b, blines in groupby(bed, key=(lambda x: x.seqid)):
        blines = list(blines)
        maxlen = max(x.end for x in blines)
        chr_lens[b] = maxlen

    for b in bed:
        accn = b.accn
        if accn == "centromere":
            centromeres[b.seqid] = b.start
        if accn in mappings:
            b.accn = mappings[accn]
        else:
            b.accn = '-'

    chr_number = len(chr_lens)
    assert chr_number == len(centromeres)

    fig = plt.figure(1, (w, h))
    root = fig.add_axes([0, 0, 1, 1])

    r = .7  # width and height of the whole chromosome set
    xstart, ystart = .15, .85
    xinterval = r / chr_number
    xwidth = xinterval * .5  # chromosome width
    max_chr_len = max(chr_lens.values())
    ratio = r / max_chr_len  # canvas / base

    # first the chromosomes
    for a, (chr, cent_position) in enumerate(sorted(centromeres.items())):
        clen = chr_lens[chr]
        xx = xstart + a * xinterval + .5 * xwidth
        yy = ystart - cent_position * ratio
        root.text(xx, ystart + .01, chr, ha="center")
        ChromosomeWithCentromere(root, xx, ystart, yy,
                ystart - clen * ratio, width=xwidth)

    chr_idxs = dict((a, i) for i, a in enumerate(sorted(chr_lens.keys())))

    alpha = .75
    # color the regions
    for chr in sorted(chr_lens.keys()):
        segment_size, excess = 0, 0
        bac_list = []
        for b in bed.sub_bed(chr):
            clen = chr_lens[chr]
            idx = chr_idxs[chr]
            klass = b.accn
            start = b.start
            end = b.end
            xx = xstart + idx * xinterval
            yystart = ystart - end * ratio
            yyend = ystart - start * ratio
            root.add_patch(Rectangle((xx, yystart), xwidth, yyend - yystart,
                fc=class_colors.get(klass, "w"), lw=0, alpha=alpha))

            if imagemap:
                """
                `segment` : size of current BAC being investigated + `excess`
                `excess`  : left-over bases from the previous BAC, as a result of
                            iterating over `winsize` regions of `segment`
                """
                if excess == 0:
                    segment_start = start
                segment = (end - start + 1) + excess
                while True:
                    if segment < winsize:
                        bac_list.append(b.accn)
                        excess = segment
                        break
                    segment_end = segment_start + winsize - 1
                    tlx, tly, brx, bry = xx, (1 - ystart) + segment_start * ratio, \
                                  xx + xwidth, (1 - ystart) + segment_end * ratio
                    print >> mapfh, '\t' + write_ImageMapLine(tlx, tly, brx, bry, \
                            w, h, dpi, chr+":"+",".join(bac_list), segment_start, segment_end)

                    segment_start += winsize
                    segment -= winsize
                    bac_list = []

        if imagemap and excess > 0:
            bac_list.append(b.accn)
            segment_end = end
            tlx, tly, brx, bry = xx, (1 - ystart) + segment_start * ratio, \
                          xx + xwidth, (1 - ystart) + segment_end * ratio
            print >> mapfh, '\t' + write_ImageMapLine(tlx, tly, brx, bry, \
                    w, h, dpi, chr+":"+",".join(bac_list), segment_start, segment_end)

    if imagemap:
        print >> mapfh, '</map>'
        mapfh.close()
        logging.debug("Image map written to `{0}`".format(mapfh.name))

    if opts.gauge:
        xstart, ystart = .9, .85
        Gauge(root, xstart, ystart - r, ystart, max_chr_len)

    # class legends, four in a row
    xstart = .1
    xinterval = .2
    xwidth = .04
    yy = .08
    for klass, cc in sorted(class_colors.items()):
        if klass == '-':
            continue
        root.add_patch(Rectangle((xstart, yy), xwidth, xwidth, fc=cc, lw=0,
            alpha=alpha))
        root.text(xstart + xwidth + .01, yy, klass, fontsize=10)
        xstart += xinterval

    empty = opts.empty
    if empty:
        root.add_patch(Rectangle((xstart, yy), xwidth, xwidth, fill=False, lw=1))
        root.text(xstart + xwidth + .01, yy, empty, fontsize=10)

    root.text(.5, .95, opts.title, fontstyle="italic", ha="center", va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    savefig(figname, dpi=dpi, iopts=iopts)


if __name__ == '__main__':
    main()
