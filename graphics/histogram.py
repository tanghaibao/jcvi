#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Use R ggplot2 library to plot histogram
"""

import os.path as op
import sys

from optparse import OptionParser

from jcvi.apps.R import RTemplate
from jcvi.apps.base import ActionDispatcher, debug
debug()

histogram_template = """
library(ggplot2)
vmin <- $vmin
vmax <- $vmax
data <- read.table('$numberfile')
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
data <- read.table('$numberfile', header=T, sep="\t")
m <- ggplot(data, aes(x=$xlabel, fill=group))
m + geom_bar(binwidth=(vmax-vmin)/$bins, position="dodge") + xlim(vmin, vmax) +
opts(title='$title')
ggsave('$outfile')
"""


def histogram(numberfile, vmin, vmax, xlabel, title, bins=50):
    """
    Generate histogram using number from numberfile, and only numbers in the
    range of (vmin, vmax)
    """
    outfile = numberfile + '.pdf'

    rtemplate = RTemplate(histogram_template, locals())
    rtemplate.run()


def histogram_multiple(numberfiles, vmin, vmax, xlabel, title, bins=50):
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
    p.add_option("--vmin", dest="vmin", default=-1e9, type="int",
            help="numbers larger than this is ignored [default: %default]")
    p.add_option("--vmax", dest="vmax", default=1e9, type="int",
            help="numbers larger than this is ignored [default: %default]")
    p.add_option("--bins", dest="bins", default=50, type="int",
            help="number of bins to plot in the histogram [default: %default]")
    p.add_option("--xlabel", dest="xlabel", default="value",
            help="label on the X-axis")
    p.add_option("--title", dest="title", default="title",
            help="title of the plot")
    opts, args = p.parse_args()

    if len(args) < 1:
        sys.exit(p.print_help())

    vmin, vmax = opts.vmin, opts.vmax
    bins = opts.bins
    xlabel, title = opts.xlabel, opts.title

    fileno = len(args)
    if fileno == 1:
        histogram(args[0], vmin, vmax, xlabel, title, bins=bins)
    else:
        histogram_multiple(args, vmin, vmax, xlabel, title, bins=bins)


if __name__ == '__main__':
    main()
