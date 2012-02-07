#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Collect gene statistics based on gff file:
Exon length, Intron length, Gene length, Exon count
"""

import sys

from optparse import OptionParser

from jcvi.formats.gff import GffLine, make_index
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('stats', 'collect gene statistics based on gff file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def stats(args):
    """
    %prog stats infile.gff

    Collect gene statistics based on gff file. There are some terminology issues
    here and so normally we call "gene" are actually mRNA, and "exon" are actually
    CDS, but they are configurable.
    """
    import GFFutils

    from jcvi.utils.cbook import SummaryStats
    from jcvi.utils.range import range_interleave

    p = OptionParser(stats.__doc__)
    p.add_option("--gene", default="mRNA",
                 help="The gene type [default: %default]")
    p.add_option("--exon", default="CDS",
                 help="The exon type [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gff_file, = args
    db_file = make_index(gff_file)
    g = GFFutils.GFFDB(db_file)
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

    print >> sys.stderr, SummaryStats(exon_lengths, title="Exon Length")
    print >> sys.stderr, SummaryStats(intron_lengths, title="Intron Length")
    print >> sys.stderr, SummaryStats(gene_lengths, title="Gene Length")
    print >> sys.stderr, SummaryStats(exon_counts, title="Exon Count")


if __name__ == '__main__':
    main()
