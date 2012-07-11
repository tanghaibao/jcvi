#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Legacy script to plot distribution of certain classes onto chromosomes. Adapted
from the script used in the Tang et al. PNAS 2010 paper, sigma figure.
"""

import sys
import os
import logging

from math import radians, ceil
from itertools import groupby
from optparse import OptionParser

import numpy as np

from jcvi.formats.bed import Bed
from jcvi.graphics.base import plt, Rectangle, Polygon, \
        CirclePolygon, _, set_image_options
from jcvi.graphics.glyph import BaseGlyph


def plot_cap(center, t, r):
    x, y = center
    return zip(x + r * np.cos(t), y + r * np.sin(t))


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


class Chromosome (object):
    def __init__(self, ax, x, y1, y2, width=.015, fc="k", fill=False, zorder=2):
        """
        Chromosome with positions given in (x, y1) => (x, y2)
        """
        pts = []
        r = width * .5
        pts += plot_cap((x, y1 - r), np.radians(range(180)), r)
        pts += [[x - r, y1 - r], [x - r, y2 + r]]
        pts += plot_cap((x, y2 + r), np.radians(range(180, 360)), r)
        pts += [[x + r, y2 + r], [x + r, y1 - r]]
        ax.add_patch(Polygon(pts, fc=fc, fill=fill, zorder=zorder))


class HorizontalChromosome (object):
    def __init__(self, ax, x1, x2, y, height=.015, ec="k", fc=None,
                       fill=False, zorder=2):
        """
        Chromosome with positions given in (x1, y) => (x2, y)
        """
        pts = self.get_pts(x1, x2, y, height)
        ax.add_patch(Polygon(pts, fill=False, ec=ec, zorder=zorder))
        if fc:
            pts = self.get_pts(x1, x2, y, height * .5)
            ax.add_patch(Polygon(pts, fc=fc, lw=0, zorder=zorder))

    def get_pts(self, x1, x2, y, height):
        r = height / (3 ** .5)

        if x2 - x1 < 2 * height:  # rectangle for small chromosomes
            return [[x1, y + r], [x1, y - r], [x2, y - r], [x2, y + r]]

        pts = []
        pts += plot_cap((x1 + r, y), np.radians(range(120, 240)), r)
        pts += [[x1 + r / 2, y - r], [x2 - r / 2, y - r]]
        pts += plot_cap((x2 - r, y), np.radians(range(-60, 60)), r)
        pts += [[x2 - r / 2, y + r], [x1 + r / 2, y + r]]

        return pts


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
    opts, args, iopts = set_image_options(p, figsize="6x6", dpi=300)

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
        mappings = dict(x.split() for x in open(mappingfile))
        classes = sorted(set(mappings.values()))
        logging.debug("A total of {0} classes found: {1}".format(len(classes),
            ','.join(classes)))
    else:
        mappings = {}
        classes = []
        logging.debug("No classes registered (no id_mappings given).")

    mycolors = "wrgbymc"
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
        root.text(xx, ystart + .01, _(chr), ha="center")
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
        tip = .008  # the ticks on the gauge bar
        extra = .006  # the offset for the unit label
        xstart, ystart = .9, .85
        yy = ystart
        gauge = int(ceil(max_chr_len / 1e6))
        mb = ratio * 1e6
        yinterval = 2 * mb
        root.plot([xstart, xstart], [yy, yy - r], 'b-', lw=2)
        for x in xrange(0, gauge, 2):
            if x % 10:
                root.plot([xstart, xstart + tip], [yy, yy], "b-")
            else:
                root.plot([xstart - tip, xstart + tip], [yy, yy], 'b-', lw=2)
                root.text(xstart + tip + extra, yy, _(x),
                        color="gray", va="center")
            yy -= yinterval
        root.text(xstart, yy - .03, _("Mb"), color="gray", va="center")

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
        root.text(xstart + xwidth + .01, yy, _(klass), fontsize=9)
        xstart += xinterval

    root.text(.5, .95, opts.title, fontstyle="italic", ha="center", va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    plt.savefig(figname, dpi=dpi)
    logging.debug("Figure saved to `{0}` {1}".format(figname, iopts))


if __name__ == '__main__':
    main()
