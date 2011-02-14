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
from optparse import OptionParser

import numpy as np

from jcvi.graphics.base import plt, ticker, Rectangle, cm, _, \
    human_size_formatter
from jcvi.formats.base import LineFile
from jcvi.formats.blast import BlastLine
from jcvi.formats.bed import Bed
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
        sizes = list(self.iter_sizes())
        ctgs, sizes = zip(*sizes)
        sizes = np.cumsum([0] + list(sizes))
        self.mapping = dict(zip(ctgs, sizes))
        self.ctgs = ctgs
        self.sizes = sizes

    def __len__(self):
        return len(self.mapping)

    @property
    def totalsize(self):
        return np.max(self.sizes)

    def iter_sizes(self):
        self.fp.seek(0)
        for row in self.fp:
            ctg, size = row.split()[:2]
            yield ctg, int(size)

    def get_position(self, ctg, pos):
        return self.mapping[ctg] + pos

    def get_breaks(self):
        for i in xrange(1, len(self)):
            yield self.ctgs[i], self.sizes[i-1], self.sizes[i]


def blastplot(blastfile, qsizes, ssizes, qbed, sbed, image_name):

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
            qi = (q.start + q.end) / 2
        else:
            qi = (b.qstart + b.qstop) / 2 

        if sorder:
            if subject not in sorder: continue
            si, s = sorder[subject]
            subject = s.seqid
            si = (s.start + s.end) / 2
        else:
            si = (b.sstart + b.sstop) / 2

        qi = qsizes.get_position(query, qi)
        si = ssizes.get_position(subject, si)
        data.append((qi, si))

    fig = plt.figure(1,(8,8))
    root = fig.add_axes([0,0,1,1]) # the whole canvas
    ax = fig.add_axes([.1,.1,.8,.8]) # the dot plot

    """
    sample_number = 5000 # only show random subset
    if len(data) > sample_number:
        data = sample(data, sample_number)
    """

    x, y = zip(*data)
    ax.scatter(x, y, s=2, lw=0)

    xsize, ysize = qsizes.totalsize, ssizes.totalsize
    logging.debug("xsize=%d ysize=%d" % (xsize, ysize))
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
    logging.debug("print image to %s" % image_name)
    plt.savefig(image_name, dpi=150)


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

    qsizes = Sizes(qsizes)
    ssizes = Sizes(ssizes)
    if qbed: qbed = Bed(qbed)
    if sbed: sbed = Bed(sbed)

    blastfile = args[0]

    image_name = op.splitext(blastfile)[0] + "." + opts.format
    blastplot(blastfile, qsizes, ssizes, qbed, sbed, image_name)

