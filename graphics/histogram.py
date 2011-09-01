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

import numpy as np

from jcvi.graphics.base import asciiplot
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

def loghistogram(bins, base=2, ascii=True, title="Counts"):
    """
    bins is a dictionary with key: log(x, base), value: counts.
    """
    x, y = [], []
    for size, number in sorted(bins.items()):
        lb, ub = base ** size, base ** (size + 1)
        x.append((lb, ub))
        y.append(number)

    asciiplot(x, y, title=title)


def get_data(filename, vmin=None, vmax=None, skip=0):
    fp = open(filename)
    for s in xrange(skip):
        fp.next()
    data = np.array([float(x) for x in fp])
    vmin = min(data) if vmin is None else vmin
    vmax = max(data) if vmax is None else vmax
    data = data[(data >= vmin) & (data <= vmax)]

    return data, vmin, vmax


def stem_leaf_plot(data, vmin, vmax, bins, digit=1, title=None):
    '''
    Generate stem and leaf plot given a collection of numbers
    '''
    assert bins > 0
    range = vmax - vmin
    if range % bins == 0:
        step = range / bins
    else:
        step = range * 1. / bins
    step = step or 1

    bins = np.arange(vmin, vmax + step, step)
    hist, bin_edges = np.histogram(data, bins=bins)
    asciiplot(bin_edges, hist, digit=digit, title=title)
    print >> sys.stderr, "Last bin ends in {0}, inclusive.".format(vmax)


def texthistogram(numberfiles, vmin, vmax, title=None, bins=20, skip=0):
    for nf in numberfiles:
        logging.debug("Import `{0}`.".format(nf))
        data, vmin, vmax = get_data(nf, vmin, vmax, skip=skip)
        stem_leaf_plot(data, vmin, vmax, bins, title=title)


def histogram(numberfile, vmin, vmax, xlabel, title,
        bins=50, skip=0, ascii=False):
    """
    Generate histogram using number from numberfile, and only numbers in the
    range of (vmin, vmax)
    """
    if ascii:
        return texthistogram([numberfile], vmin, vmax, title=title,
                bins=bins, skip=skip)

    outfile = numberfile + '.pdf'
    data, vmin, vmax = get_data(numberfile, vmin, vmax, skip=skip)

    rtemplate = RTemplate(histogram_template, locals())
    rtemplate.run()


def histogram_multiple(numberfiles, vmin, vmax, xlabel, title,
        bins=20, skip=0, ascii=False):
    """
    Generate histogram using number from numberfile, and only numbers in the
    range of (vmin, vmax). First combining multiple files.
    """
    if ascii:
        return texthistogram(numberfiles, vmin, vmax, title=title,
                bins=bins, skip=skip)

    newfile = "_".join(op.basename(x).split(".")[0] for x in numberfiles)
    suffix = op.basename(numberfiles[0]).split(".")[-1]
    newfile += "." + suffix

    fw = open(newfile, "w")
    print >> fw, "{0}\tgroup".format(xlabel)
    for f in numberfiles:
        data, va, vb = get_data(f, vmin, vmax, skip=skip)
        vmin = min(vmin, va)
        vmax = max(vmax, vb)

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
    p.add_option("--vmin", dest="vmin", default=None, type="int",
            help="minimum value, inclusive [default: %default]")
    p.add_option("--vmax", dest="vmax", default=None, type="int",
            help="maximum value, inclusive [default: %default]")
    p.add_option("--bins", dest="bins", default=20, type="int",
            help="number of bins to plot in the histogram [default: %default]")
    p.add_option("--xlabel", dest="xlabel", default="value",
            help="label on the X-axis")
    p.add_option("--title", help="title of the plot")
    p.add_option("--ascii", default=False, action="store_true",
        help="print ASCII text stem-leaf plot [default: %default]")
    opts, args = p.parse_args()

    if len(args) < 1:
        sys.exit(not p.print_help())

    skip = opts.skip
    vmin, vmax = opts.vmin, opts.vmax
    bins = opts.bins
    xlabel, title = opts.xlabel, opts.title
    title = title or args[0]

    fileno = len(args)
    if fileno == 1:
        histogram(args[0], vmin, vmax, xlabel, title,
                bins=bins, skip=skip, ascii=opts.ascii)
    else:
        histogram_multiple(args, vmin, vmax, xlabel, title,
                bins=bins, skip=skip, ascii=opts.ascii)


if __name__ == '__main__':
    main()
