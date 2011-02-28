#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog blastfile --qsizes query.sizes --ssizes subject.sizes

visualize the blastfile in a dotplot.
"""

import os.path as op
import sys
import logging
from random import sample
from optparse import OptionParser

import numpy as np

from jcvi.graphics.base import plt, ticker, Rectangle, cm, _, \
    human_size_formatter
from jcvi.formats.blast import BlastLine
from jcvi.formats.sizes import Sizes
from jcvi.formats.bed import Bed
from jcvi.apps.base import debug
debug()


def blastplot(blastfile, qsizes, ssizes, qbed, sbed, image_name, 
        lines=False, proportional=False, sampleN=5000):

    fp = open(blastfile)

    qorder = qbed.order if qbed else None
    sorder = sbed.order if sbed else None

    data = []

    for row in fp:
        b = BlastLine(row)
        query, subject = b.query, b.subject

        if qorder:
            if query not in qorder: continue
            qi, q = qorder[query]
            query = q.seqid
            qstart, qend = q.start, q.end
        else:
            qstart, qend = b.qstart, b.qstop

        if sorder:
            if subject not in sorder: continue
            si, s = sorder[subject]
            subject = s.seqid
            sstart, send = s.start, s.end
        else:
            sstart, send = b.sstart, b.sstop

        qi = qsizes.get_position(query, qstart)
        qj = qsizes.get_position(query, qend)
        si = ssizes.get_position(subject, sstart)
        sj = ssizes.get_position(subject, send)

        if None in (qi, si): continue
        data.append(((qi, qj), (si, sj)))

    if sampleN:
        if len(data) > sampleN:
            data = sample(dataN)

    if not data:
        return logging.error("no blast data imported")

    xsize, ysize = qsizes.totalsize, ssizes.totalsize
    logging.debug("xsize=%d ysize=%d" % (xsize, ysize))

    # Fix the width
    ratio = 1
    if proportional: ratio = ysize * 1./ xsize
    width = 8
    height = width * ratio
    fig = plt.figure(1, (width, height))
    root = fig.add_axes([0,0,1,1]) # the whole canvas
    ax = fig.add_axes([.1,.1,.8,.8]) # the dot plot

    if lines:
        for a, b in data:
            ax.plot(a, b, 'ro-', mfc="w", mec="r", ms=3)
    else:
        data = [(x[0], y[0]) for x, y in data]
        x, y = zip(*data)
        ax.scatter(x, y, s=2, lw=0)

    xlim = (0, xsize)
    ylim = (ysize, 0) # invert the y-axis

    xchr_labels, ychr_labels = [], []
    ignore = True # tag to mark whether to plot chromosome name (skip small ones)
    ignore_size_x = xsize * .005 
    ignore_size_y = ysize * .005

    # plot the chromosome breaks
    logging.debug("adding query breaks (%d)" % len(qsizes))
    for (seqid, beg, end) in qsizes.get_breaks():
        ignore = abs(end-beg) < ignore_size_x
        if ignore: continue
        seqid = seqid.split("_")[-1]
        try: 
            seqid = int(seqid)
            seqid = "c%d" % seqid
        except: 
            pass

        xchr_labels.append((seqid, (beg + end)/2, ignore))
        ax.plot([beg, beg], ylim, "g-", lw=1)

    logging.debug("adding subject breaks (%d)" % len(ssizes))
    for (seqid, beg, end) in ssizes.get_breaks():
        ignore = abs(end-beg) < ignore_size_y
        if ignore: continue
        seqid = seqid.split("_")[-1]
        try: 
            seqid = int(seqid)
            seqid = "c%d" % seqid
        except:
            pass

        ychr_labels.append((seqid, (beg + end)/2, ignore))
        ax.plot(xlim, [beg, beg], "g-", lw=1)

    # plot the chromosome labels
    for label, pos, ignore in xchr_labels:
        pos = .1 + pos * .8/ xsize
        if not ignore:
            root.text(pos, .91, _(label), color="b",
                va="bottom", rotation=45)

    # remember y labels are inverted
    for label, pos, ignore in ychr_labels:
        pos = .9 - pos * .8/ ysize
        if not ignore:
            root.text(.91, pos, _(label), color="b",
                ha="left", va="center")

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    to_ax_label = lambda fname: _(op.basename(fname).split(".")[0])

    # add genome names
    ax.set_xlabel(to_ax_label(qsizes.filename), size=15)
    ax.set_ylabel(to_ax_label(ssizes.filename), size=15)

    # beautify the numeric axis
    for tick in ax.get_xticklines() + ax.get_yticklines():
        tick.set_visible(False) 

    formatter = human_size_formatter
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

    plt.setp(ax.get_xticklabels() + ax.get_yticklabels(), color='gray', size=10)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    logging.debug("print image to `{0}`".format(image_name))
    plt.savefig(image_name, dpi=150)


if __name__ == "__main__":

    p = OptionParser(__doc__)
    p.add_option("--qsizes", dest="qsizes", help="path to two column qsizes file")
    p.add_option("--ssizes", dest="ssizes", help="path to two column ssizes file")
    p.add_option("--qbed", dest="qbed", help="path to qbed")
    p.add_option("--sbed", dest="sbed", help="path to sbed")
    p.add_option("--lines", dest="lines", default=False, action="store_true",
            help="plot lines for anchors instead of points [default: points]")
    p.add_option("--proportional", dest="proportional",
            default=False, action="store_true", 
            help="make the image width/height equal to seqlen ratio")
    p.add_option("--format", dest="format", default="png",
            help="generate image of format (png, pdf, ps, eps, svg, etc.)"
            "[default: %default]")

    opts, args = p.parse_args()

    qsizes, ssizes = opts.qsizes, opts.ssizes
    qbed, sbed = opts.qbed, opts.sbed
    lines = opts.lines
    proportional = opts.proportional

    if not (len(args) == 1 and qsizes and ssizes):
        sys.exit(p.print_help())

    qsizes = Sizes(qsizes)
    ssizes = Sizes(ssizes)
    if qbed: qbed = Bed(qbed)
    if sbed: sbed = Bed(sbed)

    blastfile = args[0]

    image_name = op.splitext(blastfile)[0] + "." + opts.format
    blastplot(blastfile, qsizes, ssizes, qbed, sbed, image_name, 
            lines=lines, proportional=proportional)

