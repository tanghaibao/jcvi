#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog blastfile --qsizes query.sizes --sbed subject.sizes

visualize the blastfile in a dotplot.
"""

import os.path as op
import sys
import logging
from random import sample
from itertools import groupby
from optparse import OptionParser

import numpy as np

from jcvi.graphics.base import plt, ticker, Rectangle, cm, _
from jcvi.formats.base import LineFile
from jcvi.formats.blast import BlastLine
from jcvi.formats.bed import Bed
from jcvi.algorithms.synteny import batch_scan
from jcvi.apps.base import debug
debug()


class Sizes (LineFile):
    """
    Two-column file,
    contigID<tab>size
    """
    def __init__(self, filename):
        super(Sizes, self).__init__(filename)
        self.fp = open(filename)

    def iter_sizes(self):
        self.fp.seek(0)
        for row in self.fp:
            ctg, size = atoms.split()[:2]
            yield ctg, int(size)


def get_breaks(bed):
    # get chromosome break positions
    simple_bed = bed.simple_bed
    for seqid, ranks in groupby(simple_bed, key=lambda x:x[0]):
        ranks = list(ranks)
        # chromosome, extent of the chromosome
        yield seqid, ranks[0][1], ranks[-1][1]


def blastplot(blastfile, qsizes, ssizes, qbed, sbed, image_name):

    fp = open(blastfile)

    qorder = qbed.order
    sorder = sbed.order

    data = []
    logging.debug("normalize the values to [%.1f, %.1f]" % (vmin, vmax))

    for row in fp:
        atoms = row.split()
        if len(atoms) < 3: continue
        query, subject, value = row.split()[:3]
        try: value = float(value)
        except: value = vmin 

        if value < vmin: value = vmin
        if value > vmax: value = vmax

        if query not in qorder: 
            #logging.warning("ignore %s" % query)
            continue
        if subject not in sorder: 
            #logging.warning("ignore %s" % subject)
            continue

        qi, q = qorder[query]
        si, s = sorder[subject]
        data.append((qi, si, vmax-value))

    fig = plt.figure(1,(8,8))
    root = fig.add_axes([0,0,1,1]) # the whole canvas
    ax = fig.add_axes([.1,.1,.8,.8]) # the dot plot

    sample_number = 5000 # only show random subset
    if len(data) > sample_number:
        data = sample(data, sample_number)

    # the data are plotted in this order, the least value are plotted
    # last for aesthetics
    data.sort(key=lambda x: -x[2])

    x, y, c = zip(*data)
    ax.scatter(x, y, c=c, s=2, lw=0, cmap=default_cm, 
            vmin=vmin, vmax=vmax)

    xsize, ysize = len(qbed), len(sbed)
    xlim = (0, xsize)
    ylim = (ysize, 0) # invert the y-axis

    xchr_labels, ychr_labels = [], []
    ignore = True # tag to mark whether to plot chromosome name (skip small ones)
    ignore_size = 100
    # plot the chromosome breaks
    for (seqid, beg, end) in get_breaks(qbed):
        ignore = abs(end-beg) < ignore_size 
        seqid = seqid.split("_")[-1]
        try: 
            seqid = int(seqid)
            seqid = "c%d" % seqid
        except: 
            pass

        xchr_labels.append((seqid, (beg + end)/2, ignore))
        ax.plot([beg, beg], ylim, "g-", lw=1)

    for (seqid, beg, end) in get_breaks(sbed):
        ignore = abs(end-beg) < ignore_size 
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

    # create a diagonal to separate mirror image for self comparison
    if is_self:
        ax.plot(xlim, ylim, 'm-', alpha=.5, lw=2)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    to_ax_label = lambda fname: _(op.basename(fname).split(".")[0])

    # add genome names
    ax.set_xlabel(to_ax_label(qbed.filename), size=15)
    ax.set_ylabel(to_ax_label(sbed.filename), size=15)

    # beautify the numeric axis
    for tick in ax.get_xticklines() + ax.get_yticklines():
        tick.set_visible(False) 

    formatter = ticker.FuncFormatter(lambda x, pos: _("%dK" % (int(x)/1000)))
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

    plt.setp(ax.get_xticklabels() + ax.get_yticklabels(), color='gray', size=10)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    logging.debug("print image to %s" % image_name)
    plt.savefig(image_name, dpi=1000)


if __name__ == "__main__":

    p = OptionParser(__doc__)
    p.add_option("--qsizes", dest="qsizes", help="path to two column qsizes file")
    p.add_option("--ssizes", dest="ssizes", help="path to two column ssizes file")
    p.add_option("--qbed", dest="qbed", help="path to qbed")
    p.add_option("--sbed", dest="sbed", help="path to sbed")
    p.add_option("--format", dest="format", default="png",
            help="generate image of format (png, pdf, ps, eps, svg, etc.)"
            "[default: %default]")

    opts, args = p.parse_args()

    qsizes, ssizes = opts.qsizes, opts.ssizes
    qbed, sbed = opts.qbed, opts.sbed
    if not (len(args) == 1 and qsizes and ssizes):
        sys.exit(p.print_help())

    qsizes = dict(Sizes(qsizes))
    ssizes = dict(Sizes(ssizes))
    if qbed: qbed = Bed(qbed)
    if sbed: sbed = Bed(sbed)

    blastfile = args[0]

    image_name = op.splitext(blastfile)[0] + "." + opts.format
    blastplot(blastfile, qsizes, ssizes, qbed, sbed, image_name)

