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
from jcvi.graphics.base import plt, Rectangle, Polygon, CirclePolygon, _

al = .5
z = 2


def plot_cap(center, t, r):
    x, y = center
    return zip(x + r * np.cos(t), y + r * np.sin(t))


class Chromosome (object):
    def __init__(self, ax, x, y1, y2, width=.015, cl="k", zorder=z):
        """
        Chromosome with positions given in (x, y1) => (x, y2)
        """
        pts = []
        r = width * .5
        pts += plot_cap((x, y1 - r), np.radians(range(180)), r)
        pts += [[x - r, y1 - r], [x - r, y2 + r]]
        pts += plot_cap((x, y2 + r), np.radians(range(180, 360)), r)
        pts += [[x + r, y2 + r], [x + r, y1 - r]]
        ax.add_patch(Polygon(pts, fc=cl, fill=False, zorder=z))


class ChromsomeWithCentromere (object):
    def __init__(self, ax, x, y1, y2, y3, width=.015, cl="k", zorder=z):
        """
        Chromosome with centromeres at y2 position
        """
        pts = []
        r = width * .5
        pts += plot_cap((x, y1 - r), np.radians(range(180)), r)
        pts += [[x - r, y1 - r], [x - r, y2 + r]]
        pts += plot_cap((x, y2 + r), np.radians(range(180, 360)), r)
        pts += [[x + r, y2 + r], [x + r, y1 - r]]
        ax.add_patch(Polygon(pts, fc=cl, fill=False, zorder=z))
        pts = []
        pts += plot_cap((x, y2 - r), np.radians(range(180)), r)
        pts += [[x - r, y2 - r], [x - r, y3 + r]]
        pts += plot_cap((x, y3 + r), np.radians(range(180, 360)), r)
        pts += [[x + r, y3 + r], [x + r, y2 - r]]
        ax.add_patch(Polygon(pts, fc=cl, fill=False, zorder=z))
        ax.add_patch(CirclePolygon((x, y2), radius=r * .5,
            fc="k", ec="k", zorder=z))


def main():
    """
    %prog bedfile id_mappings

    Takes a bedfile that contains the coordinates of features to plot on the
    chromosomes, and `id_mappings` file that map the ids to certain class. Each
    class will get assigned a unique color.
    """
    p = OptionParser(main.__doc__)
    p.add_option("--gauge", dest="gauge", default=False, action="store_true",
            help="draw a gauge with size label [default: %default]")
    opts, args = p.parse_args()

    if len(args) != 2:
        sys.exit(p.print_help())

    bedfile, mappingfile = args
    prefix = bedfile.rsplit(".", 1)[0]
    figname = prefix + ".pdf"

    mappings = dict(x.split() for x in open(mappingfile))
    classes = sorted(set(mappings.values()))
    logging.debug("A total of {0} classes found: {1}".format(len(classes),
        ','.join(classes)))

    #mycolors = "A6C953|2F7EC1|019D75|00B2EC|EA2790|F6AE26|EA2D2D".split("|")
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

    fig = plt.figure(1, (6, 6))
    root = fig.add_axes([0, 0, 1, 1])

    r = .7  # width and height of the whole chromosome set
    xstart, ystart = .15, .15
    xinterval = r / chr_number
    xwidth = xinterval * .5  # chromosome width
    max_chr_len = max(chr_lens.values())
    ratio = r / max_chr_len  # canvas / base

    # first the chromosomes
    for a, (chr, cent_position) in enumerate(sorted(centromeres.items())):
        clen = chr_lens[chr]
        xx = xstart + a * xinterval + .5 * xwidth
        yy = ystart + (clen - cent_position) * ratio
        root.text(xx, ystart - .01, _(chr), va="top", ha="center")
        ChromsomeWithCentromere(root, xx, ystart + clen * ratio, yy,
                ystart, width=xwidth, zorder=z)

    chr_idxs = dict((a, i) for i, a in enumerate(sorted(chr_lens.keys())))

    alpha = .75
    # color the regions
    for b in bed:
        chr = b.seqid
        clen = chr_lens[chr]
        idx = chr_idxs[chr]
        xx = xstart + idx * xinterval
        start = b.start
        end = b.end
        klass = b.accn
        yystart = ystart + (clen - end) * ratio
        yyend = ystart + (clen - start) * ratio
        root.add_patch(Rectangle((xx, yystart), xwidth, yyend - yystart,
            fc=class_colors[klass], lw=0, alpha=alpha))

    if opts.gauge:
        tip = .008  # the ticks on the gauge bar
        extra = .006  # the offset for the unit label
        xstart, ystart = .9, .15
        yy = ystart
        gauge = int(ceil(max_chr_len / 1e6))
        mb = ratio * 1e6
        yinterval = 2 * mb
        root.plot([xstart, xstart], [yy, yy + r], 'b-', lw=2)
        for x in xrange(0, gauge, 2):
            if x % 10:
                root.plot([xstart, xstart + tip], [yy, yy], "b-")
            else:
                root.plot([xstart - tip, xstart + tip], [yy, yy], 'b-', lw=2)
                root.text(xstart + tip + extra, yy, _(x),
                        color="gray", va="center")
            yy += yinterval
        root.text(xstart, ystart - .03, _("Mb"), color="gray", va="center")

    # class legends, four in a row
    xstart = .1
    xinterval = .2
    xwidth = .04
    yy = .9
    for klass, cc in sorted(class_colors.items()):
        if klass == '-':
            continue
        root.add_patch(Rectangle((xstart, yy), xwidth, xwidth, fc=cc, lw=0,
            alpha=alpha))
        root.text(xstart + xwidth + .01, yy, _(klass), fontsize=9)
        xstart += xinterval

    label = "Medicago truncatula v3.5"
    root.text(.5, .05, label, fontstyle="italic", ha="center", va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    plt.savefig(figname, dpi=300)
    logging.debug("Figure saved to `{0}`".format(figname))


if __name__ == '__main__':
    main()
