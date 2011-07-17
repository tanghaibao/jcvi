#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Use R ggplot2 library to plot histogram, also contains an ASCII histogram (use
--text) when invoking histogram().
"""

import os.path as op
import sys
import logging

from optparse import OptionParser

from jcvi.apps.R import RTemplate
from jcvi.apps.base import ActionDispatcher, debug
debug()

histogram_template = """
library(ggplot2)
vmin <- $vmin
vmax <- $vmax
data <- read.table('$numberfile', skip=$skip)
data <- data[data >= vmin]
data <- data[data <= vmax]
data <- data.frame($xlabel=data)
qplot($xlabel, data=data, geom='histogram', main='$title', binwidth=(vmax-vmin)/$bins)
ggsave('$outfile')
"""

histogram_multiple_template = """
library(ggplot2)
vmin <- $vmin
vmax <- $vmax
data <- read.table('$numberfile', header=T, sep="\t", skip=$skip)
m <- ggplot(data, aes(x=$xlabel, fill=group))
m + geom_bar(binwidth=(vmax-vmin)/$bins, position="dodge") + xlim(vmin, vmax) +
opts(title='$title')
ggsave('$outfile')
"""

def stem_leaf_plot(data, bins, char="="):
    '''
    Generate stem and leaf plot given a collection of numbers
    '''
    import numpy as np

    assert bins > 0
    ma, mb = min(data), max(data)
    step = ((mb - ma) / bins) or 1

    bins = np.arange(ma, mb + step, step)

    hist, bin_edges = np.histogram(data, bins=bins)
    width = 50  # the textwidth (height) of the distribution
    hist = hist * width / hist.sum()
    for b, h in zip(bin_edges, hist):
        pct = "{0:.1f}".format(b)
        print >> sys.stderr, "{0}|{1}".format(pct.rjust(10), char * h)


def texthistogram(numberfiles, vmin, vmax, bins=50, skip=0):
    for nf in numberfiles:
        fp = open(nf)
        logging.debug("Import `{0}`.".format(nf))
        for s in xrange(skip):
            fp.next()
        data = [float(x) for x in fp]
        stem_leaf_plot(data, bins)


def histogram(numberfile, vmin, vmax, xlabel, title, bins=50, skip=0):
    """
    Generate histogram using number from numberfile, and only numbers in the
    range of (vmin, vmax)
    """
    outfile = numberfile + '.pdf'

    rtemplate = RTemplate(histogram_template, locals())
    rtemplate.run()


def histogram_multiple(numberfiles, vmin, vmax, xlabel, title, bins=50, skip=0):
    """
    Generate histogram using number from numberfile, and only numbers in the
    range of (vmin, vmax). First combining multiple files.
    """
    newfile = "_".join(op.basename(x).split(".")[0] for x in numberfiles)
    suffix = op.basename(numberfiles[0]).split(".")[-1]
    newfile += "." + suffix

    fw = open(newfile, "w")
    print >> fw, "{0}\tgroup".format(xlabel)
    for f in numberfiles:
        fp = open(f)
        tag = op.basename(f).split(".")[0]
        for row in fp:
            val = row.strip()
            print >> fw, "\t".join((val, tag))
    fw.close()

    numberfile = newfile
    outfile = numberfile + '.pdf'
    rtemplate = RTemplate(histogram_multiple_template, locals())
    rtemplate.run()


def main():
    """
    %prog numbers1.txt number2.txt ...

    Print histogram of the data files. The data files contain one number per
    line. If more than one file is inputted, the program will combine the
    histograms into the same plot.
    """
    p = OptionParser(main.__doc__)
    p.add_option("--skip", default=0, type="int",
            help="skip the first several lines [default: %default]")
    p.add_option("--vmin", dest="vmin", default=0, type="int",
            help="numbers larger than this is ignored [default: %default]")
    p.add_option("--vmax", dest="vmax", default=100000, type="int",
            help="numbers larger than this is ignored [default: %default]")
    p.add_option("--bins", dest="bins", default=50, type="int",
            help="number of bins to plot in the histogram [default: %default]")
    p.add_option("--xlabel", dest="xlabel", default="value",
            help="label on the X-axis")
    p.add_option("--title", help="title of the plot")
    p.add_option("--text", default=False, action="store_true",
        help="print ASCII text stem-leaf plot [default: %default]")
    opts, args = p.parse_args()

    if len(args) < 1:
        sys.exit(not p.print_help())

    skip = opts.skip
    vmin, vmax = opts.vmin, opts.vmax
    bins = opts.bins
    xlabel, title = opts.xlabel, opts.title
    title = title or args[0]

    if opts.text:
        texthistogram(args, vmin, vmax, bins=bins, skip=skip)
        return

    fileno = len(args)
    if fileno == 1:
        histogram(args[0], vmin, vmax, xlabel, title, bins=bins, skip=skip)
    else:
        histogram_multiple(args, vmin, vmax, xlabel, title, bins=bins, skip=skip)


if __name__ == '__main__':
    main()
