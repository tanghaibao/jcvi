#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Collect gene statistics based on gff file:
Exon length, Intron length, Gene length, Exon count
"""

import os.path as op
import sys
import logging

from jcvi.utils.cbook import SummaryStats, percentage
from jcvi.utils.table import tabulate
from jcvi.formats.gff import make_index
from jcvi.formats.base import DictFile
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, need_update


metrics = ("Exon_Length", "Intron_Length", "Gene_Length", "Exon_Count")


class GeneStats (object):

    def __init__(self, feat, conf_class, transcript_sizes, exons):
        self.fid = feat.id
        self.conf_class = conf_class
        self.num_exons = len(exons)
        self.num_transcripts = len(transcript_sizes)
        self.locus_size = feat.stop - feat.start + 1
        self.cum_transcript_size = sum(transcript_sizes)
        self.cum_exon_size = sum((stop - start + 1) \
                        for (c, start, stop) in exons)

    def __str__(self):
        return "\t".join(str(x) for x in (self.fid, self.conf_class,
                         self.num_exons, self.num_transcripts,
                         self.locus_size,
                         self.cum_transcript_size,
                         self.cum_exon_size))


def main():

    actions = (
        ('stats', 'collect gene statistics based on gff file'),
        ('genestats', 'print detailed gene statistics'),
        ('summary', 'print gene statistics table'),
        ('histogram', 'plot gene statistics based on output of stats'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def genestats(args):
    """
    %prog genestats gffile

    Print summary stats, including:
    - Number of genes
    - Number of single-exon genes
    - Number of multi-exon genes
    - Number of distinct exons
    - Number of genes with alternative transcript variants
    - Number of predicted transcripts
    - Mean number of distinct exons per gene
    - Mean number of transcripts per gene
    - Mean gene locus size (first to last exon)
    - Mean transcript size (UTR, CDS)
    - Mean exon size

    Stats modeled after barley genome paper Table 1.
    A physical, genetic and functional sequence assembly of the barley genome
    """
    p = OptionParser(genestats.__doc__)
    p.add_option("--groupby", default="conf_class",
                 help="Print separate stats groupby")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gff_file, = args
    gb = opts.groupby
    g = make_index(gff_file)

    tf = "transcript.sizes"
    if need_update(gff_file, tf):
        fw = open(tf, "w")
        for feat in g.features_of_type("mRNA"):
            fid = feat.id
            conf_class = feat.attributes.get(gb, "all")
            tsize = sum((c.stop - c.start + 1) for c in g.children(fid, 1) \
                             if c.featuretype == "exon")
            print >> fw, "\t".join((fid, str(tsize), conf_class))
        fw.close()

    tsizes = DictFile(tf, cast=int)
    conf_classes = DictFile(tf, valuepos=2)
    logging.debug("A total of {0} transcripts populated.".format(len(tsizes)))

    genes = []
    for feat in g.features_of_type("gene"):
        fid = feat.id
        transcripts = [c.id for c in g.children(fid, 1) \
                         if c.featuretype == "mRNA"]
        transcript_sizes = [tsizes[x] for x in transcripts]
        exons = set((c.chrom, c.start, c.stop) for c in g.children(fid, 2) \
                         if c.featuretype == "exon")
        conf_class = conf_classes[transcripts[0]]
        gs = GeneStats(feat, conf_class, transcript_sizes, exons)
        genes.append(gs)

    r = {}  # Report
    distinct_groups = set(conf_classes.values())
    for g in distinct_groups:
        num_genes = num_single_exon_genes = num_multi_exon_genes = 0
        num_genes_with_alts = num_transcripts = num_exons = max_transcripts = 0
        cum_locus_size = cum_transcript_size = cum_exon_size = 0
        for gs in genes:
            if gs.conf_class != g:
                continue
            num_genes += 1
            if gs.num_exons == 1:
                num_single_exon_genes += 1
            else:
                num_multi_exon_genes += 1
            num_exons += gs.num_exons
            if gs.num_transcripts > 1:
                num_genes_with_alts += 1
            if gs.num_transcripts > max_transcripts:
                max_transcripts = gs.num_transcripts
            num_transcripts += gs.num_transcripts
            cum_locus_size += gs.locus_size
            cum_transcript_size += gs.cum_transcript_size
            cum_exon_size += gs.cum_exon_size

        mean_num_exons = num_exons * 1. / num_genes
        mean_num_transcripts = num_transcripts * 1. / num_genes
        mean_locus_size = cum_locus_size * 1. / num_genes
        mean_transcript_size = cum_transcript_size * 1. / num_transcripts
        mean_exon_size = cum_exon_size * 1. / num_exons

        r[("Number of genes", g)] = num_genes
        r[("Number of single-exon genes", g)] = \
            percentage(num_single_exon_genes, num_genes, mode=1)
        r[("Number of multi-exon genes", g)] = \
            percentage(num_multi_exon_genes, num_genes, mode=1)
        r[("Number of distinct exons", g)] = num_exons
        r[("Number of genes with alternative transcript variants", g)] = \
            percentage(num_genes_with_alts, num_genes, mode=1)
        r[("Number of predicted transcripts", g)] = num_transcripts
        r[("Mean number of distinct exons per gene", g)] = mean_num_exons
        r[("Mean number of transcripts per gene", g)] = mean_num_transcripts
        r[("Max number of transcripts per gene", g)] = max_transcripts
        r[("Mean gene locus size (first to last exon)", g)] = mean_locus_size
        r[("Mean transcript size (UTR, CDS)", g)] = mean_transcript_size
        r[("Mean exon size", g)] = mean_exon_size

    print >> sys.stderr, tabulate(r)


def summary(args):
    """
    %prog summary *.gff

    Print gene statistics table.
    """
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

    With data written to disk then you can run %prog histogram
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
