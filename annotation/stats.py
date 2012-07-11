#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Collect gene statistics based on gff file:
Exon length, Intron length, Gene length, Exon count
"""

import os.path as op
import sys
import logging

from optparse import OptionParser

from jcvi.utils.cbook import SummaryStats
from jcvi.formats.gff import GffLine, make_index
from jcvi.apps.base import ActionDispatcher, debug, mkdir
debug()


metrics = ("Exon_Length", "Intron_Length", "Gene_Length", "Exon_Count")


def main():

    actions = (
        ('stats', 'collect gene statistics based on gff file'),
        ('histogram', 'plot gene statistics based on output of stats'),
        ('summary', 'print gene statistics table'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def summary(args):
    """
    %prog summary *.gff

    Print gene statistics table.
    """
    from jcvi.utils.table import tabulate

    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    gff_files = args
    for metric in metrics:
        logging.debug("Parsing files in `{0}`..".format(metric))

        table = {}
        for x in gff_files:
            pf = op.basename(x).split(".")[0]
            numberfile = op.join(metric, pf + ".txt")
            ar = [int(x.strip()) for x in open(numberfile)]
            sum = SummaryStats(ar).todict().items()
            keys, vals = zip(*sum)
            keys = [(pf, x) for x in keys]
            table.update(dict(zip(keys, vals)))

        print >> sys.stderr, tabulate(table)


def histogram(args):
    """
    %prog histogram *.gff

    Plot gene statistics based on output of stats. For each gff file, look to
    see if the metrics folder (i.e. Exon_Length) contains the data and plot
    them.
    """
    from jcvi.graphics.histogram import histogram_multiple

    p = OptionParser(histogram.__doc__)
    p.add_option("--bins", dest="bins", default=40, type="int",
            help="number of bins to plot in the histogram [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    gff_files = args
    # metrics = ("Exon_Length", "Intron_Length", "Gene_Length", "Exon_Count")
    colors = ("red", "green", "blue", "black")
    vmaxes = (1000, 1000, 4000, 20)
    xlabels = ("bp", "bp", "bp", "number")
    for metric, color, vmax, xlabel in zip(metrics, colors, vmaxes, xlabels):
        logging.debug("Parsing files in `{0}`..".format(metric))
        numberfiles = [op.join(metric, op.basename(x).split(".")[0] + ".txt") \
                        for x in gff_files]

        histogram_multiple(numberfiles, 0, vmax, xlabel, metric,
                       bins=opts.bins, facet=True, fill=color,
                       prefix=metric + ".")


def stats(args):
    """
    %prog stats infile.gff

    Collect gene statistics based on gff file. There are some terminology issues
    here and so normally we call "gene" are actually mRNA, and sometimes "exon"
    are actually CDS, but they are configurable.

    Use --txt to send the numbers to text file in four separate folders,
    corresponding to the four metrics:

    Exon length, Intron length, Gene length, Exon count

    With data written to disk then you can run $prog histogram
    """
    from jcvi.utils.range import range_interleave

    p = OptionParser(stats.__doc__)
    p.add_option("--gene", default="mRNA",
                 help="The gene type [default: %default]")
    p.add_option("--exon", default="CDS",
                 help="The exon type [default: %default]")
    p.add_option("--txt", default=False, action="store_true",
                 help="Print out numbers for further analyses [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gff_file, = args
    g = make_index(gff_file)
    exon_lengths = []
    intron_lengths = []
    gene_lengths = []
    exon_counts = []
    for feat in g.features_of_type(opts.gene):
        exons = []
        for c in g.children(feat.id, 1):
            if c.featuretype != opts.exon:
                continue
            exons.append((c.chrom, c.start, c.stop))
        introns = range_interleave(exons)
        feat_exon_lengths = [(stop - start + 1) for (chrom, start, stop) in exons]
        feat_intron_lengths = [(stop - start + 1) for (chrom, start, stop) in introns]
        exon_lengths += feat_exon_lengths
        intron_lengths += feat_intron_lengths
        gene_lengths.append(sum(feat_exon_lengths))
        exon_counts.append(len(feat_exon_lengths))

    a = SummaryStats(exon_lengths)
    b = SummaryStats(intron_lengths)
    c = SummaryStats(gene_lengths)
    d = SummaryStats(exon_counts)
    for x, title in zip((a, b, c, d), metrics):
        x.title = title
        print >> sys.stderr, x

    if not opts.txt:
        return

    prefix = gff_file.split(".")[0]
    for x in (a, b, c, d):
        dirname = x.title
        mkdir(dirname)
        txtfile = op.join(dirname, prefix + ".txt")
        x.tofile(txtfile)


if __name__ == '__main__':
    main()
