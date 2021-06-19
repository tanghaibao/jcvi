#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Use R ggplot2 library to plot histogram, also contains an ASCII histogram (use
--text) when invoking histogram().
"""
import os.path as op
import sys
import logging

import numpy as np

from math import log, ceil
from collections import defaultdict

from jcvi.formats.base import DictFile
from jcvi.graphics.base import asciiplot, quickplot
from jcvi.apps.r import RTemplate
from jcvi.apps.base import OptionParser


histogram_header = """
library(ggplot2)
vmin <- $vmin
vmax <- $vmax
data <- read.table('$numberfile', skip=$skip)
data <- data[data >= vmin]
data <- data[data <= vmax]
data <- data.frame($xlabel=data)
m <- ggplot(data, aes(x=$xlabel)) +
     theme(plot.title=element_text(size=11, colour="darkblue"))
"""

histogram_template = (
    histogram_header
    + """
m + geom_histogram(colour="darkgreen", fill="$fill", binwidth=(vmax-vmin)/$bins) +
labs(title='$title')
ggsave('$outfile')
"""
)

histogram_log_template = (
    histogram_header
    + """
library(scales)
m + geom_histogram(colour="darkgreen", fill="$fill", binwidth=0.33) +
labs(title='$title') +
scale_x_continuous(trans=log${base}_trans())
ggsave('$outfile')
"""
)

histogram_multiple_template = """
library(ggplot2)
vmin <- $vmin
vmax <- $vmax
data <- read.table('$numberfile', header=T, sep="\t", skip=$skip)
"""

histogram_multiple_template_a = (
    histogram_multiple_template
    + """
m <- ggplot(data, aes(x=$xlabel, fill=grp))
m + geom_bar(binwidth=(vmax-vmin)/$bins, position="dodge") +
xlim(vmin, vmax) +
labs(title='$title') +
ggsave('$outfile')
"""
)

histogram_multiple_template_b = (
    histogram_multiple_template
    + """
m <- ggplot(data, aes(x=$xlabel))
m + geom_histogram(colour="darkgreen", fill="$fill", binwidth=(vmax-vmin)/$bins) +
xlim(vmin, vmax) +
labs(title='$title') +
facet_wrap(~grp)
ggsave('$outfile')
"""
)


def loghistogram(data, base=2, title="Counts", summary=False):
    """
    bins is a dictionary with key: log(x, base), value: counts.
    """
    from jcvi.utils.cbook import percentage

    if summary:
        unique = len(data)
        total = sum(data)

        # Print out a distribution
        print("Unique: {0}".format(percentage(unique, total)), file=sys.stderr)

    bins = defaultdict(int)
    for d in data:
        logd = int(log(d, base))
        bins[logd] += 1

    x, y = [], []
    for size, number in sorted(bins.items()):
        lb, ub = base ** size, base ** (size + 1)
        x.append((lb, ub))
        y.append(number)

    asciiplot(x, y, title=title)


def get_data(filename, vmin=None, vmax=None, skip=0, col=0):
    from jcvi.utils.cbook import SummaryStats

    fp = open(filename)
    # Determine the data type
    for s in range(skip):
        next(fp)
    for row in fp:
        ntype = float if "." in row else int
        break

    fp = open(filename)
    for s in range(skip):
        next(fp)

    data = np.array([ntype(x.split()[col]) for x in fp])
    s = SummaryStats(data, title=filename)
    print(s, file=sys.stderr)

    vmin = min(data) if vmin is None else vmin
    vmax = max(data) if vmax is None else vmax
    data = data[(data >= vmin) & (data <= vmax)]

    return data, vmin, vmax


def stem_leaf_plot(data, vmin, vmax, bins, digit=1, title=None):
    """
    Generate stem and leaf plot given a collection of numbers
    """
    assert bins > 0
    range = vmax - vmin
    step = range * 1.0 / bins
    if isinstance(range, int):
        step = int(ceil(step))

    step = step or 1

    bins = np.arange(vmin, vmax + step, step)
    hist, bin_edges = np.histogram(data, bins=bins)
    # By default, len(bin_edges) = len(hist) + 1
    bin_edges = bin_edges[: len(hist)]
    asciiplot(bin_edges, hist, digit=digit, title=title)
    print("Last bin ends in {0}, inclusive.".format(vmax), file=sys.stderr)

    return bin_edges, hist


def texthistogram(numberfiles, vmin, vmax, title=None, bins=20, skip=0, col=0, base=0):

    for nf in numberfiles:
        logging.debug("Import `{0}`.".format(nf))
        data, vmin, vmax = get_data(nf, vmin, vmax, skip=skip, col=col)
        if base:
            loghistogram(data, base=base, title=title)
        else:
            stem_leaf_plot(data, vmin, vmax, bins, title=title)


def histogram(
    numberfile,
    vmin,
    vmax,
    xlabel,
    title,
    outfmt="pdf",
    bins=50,
    skip=0,
    col=0,
    ascii=False,
    base=0,
    fill="white",
):
    """
    Generate histogram using number from numberfile, and only numbers in the
    range of (vmin, vmax)
    """
    if ascii:
        return texthistogram(
            [numberfile],
            vmin,
            vmax,
            title=title,
            bins=bins,
            skip=skip,
            col=col,
            base=base,
        )

    data, vmin, vmax = get_data(numberfile, vmin, vmax, skip=skip, col=col)
    outfile = (
        numberfile + ".base{0}.{1}".format(base, outfmt)
        if base
        else numberfile + ".pdf"
    )
    template = histogram_log_template if base else histogram_template
    rtemplate = RTemplate(template, locals())
    rtemplate.run()


def histogram_multiple(
    numberfiles,
    vmin,
    vmax,
    xlabel,
    title,
    outfmt="pdf",
    tags=None,
    bins=20,
    skip=0,
    ascii=False,
    facet=False,
    fill="white",
    prefix="",
):
    """
    Generate histogram using number from numberfile, and only numbers in the
    range of (vmin, vmax). First combining multiple files.
    """
    if ascii:
        return texthistogram(numberfiles, vmin, vmax, title=title, bins=bins, skip=skip)

    newfile = "_".join(op.basename(x).split(".")[0] for x in numberfiles)

    fw = open(newfile, "w")
    print("{0}\tgrp".format(xlabel), file=fw)

    if tags:
        tags = tags.split(",")

    for i, f in enumerate(numberfiles):
        data, va, vb = get_data(f, vmin, vmax, skip=skip)
        vmin = min(vmin, va)
        vmax = max(vmax, vb)

        fp = open(f)
        if tags:
            tag = tags[i]
        else:
            tag = op.basename(f).rsplit(".", 1)[0]
        for row in fp:
            val = row.strip()
            print("\t".join((val, tag)), file=fw)
    fw.close()

    numberfile = newfile
    outfile = numberfile + "." + outfmt
    if prefix:
        outfile = prefix + outfile
    htemplate = (
        histogram_multiple_template_b if facet else histogram_multiple_template_a
    )
    rtemplate = RTemplate(htemplate, locals())
    rtemplate.run()


def main():
    """
    %prog numbers1.txt number2.txt ...

    Print histogram of the data files. The data files contain one number per
    line. If more than one file is inputted, the program will combine the
    histograms into the same plot.
    """
    allowed_format = ("emf", "eps", "pdf", "png", "ps", "raw", "rgba", "svg", "svgz")
    p = OptionParser(main.__doc__)
    p.add_option("--skip", default=0, type="int", help="skip the first several lines")
    p.add_option("--col", default=0, type="int", help="Get the n-th column")
    p.set_histogram()
    p.add_option(
        "--tags",
        dest="tags",
        default=None,
        help="tags for data if multiple input files, comma sep",
    )
    p.add_option(
        "--ascii",
        default=False,
        action="store_true",
        help="print ASCII text stem-leaf plot",
    )
    p.add_option(
        "--base",
        default="0",
        choices=("0", "2", "10"),
        help="use logarithm axis with base, 0 to disable",
    )
    p.add_option(
        "--facet",
        default=False,
        action="store_true",
        help="place multiple histograms side-by-side",
    )
    p.add_option("--fill", default="white", help="color of the bin")
    p.add_option(
        "--format",
        default="pdf",
        choices=allowed_format,
        help="Generate image of format",
    )
    p.add_option(
        "--quick",
        default=False,
        action="store_true",
        help="Use quick plot, assuming bins are already counted",
    )
    p.add_option(
        "--noprintstats",
        default=False,
        action="store_true",
        help="Write basic stats when using --quick",
    )
    opts, args = p.parse_args()

    if len(args) < 1:
        sys.exit(not p.print_help())

    skip = opts.skip
    vmin, vmax = opts.vmin, opts.vmax
    bins = opts.bins
    xlabel, title = opts.xlabel, opts.title
    title = title or args[0]
    base = int(opts.base)
    fileno = len(args)

    if opts.quick:
        assert fileno == 1, "Single input file expected using --quick"
        filename = args[0]
        figname = filename.rsplit(".", 1)[0] + ".pdf"
        data = DictFile(filename, keycast=int, cast=int)
        quickplot(
            data,
            vmin,
            vmax,
            xlabel,
            title,
            figname=figname,
            print_stats=(not opts.noprintstats),
        )
        return

    if fileno == 1:
        histogram(
            args[0],
            vmin,
            vmax,
            xlabel,
            title,
            outfmt=opts.format,
            bins=bins,
            skip=skip,
            ascii=opts.ascii,
            base=base,
            fill=opts.fill,
            col=opts.col,
        )
    else:
        histogram_multiple(
            args,
            vmin,
            vmax,
            xlabel,
            title,
            outfmt=opts.format,
            tags=opts.tags,
            bins=bins,
            skip=skip,
            ascii=opts.ascii,
            facet=opts.facet,
            fill=opts.fill,
        )


if __name__ == "__main__":
    main()
