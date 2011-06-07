#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Benchmark de novo assemblies by mapping to known sequence set.
"""

import sys

from optparse import OptionParser

from jcvi.algorithms.supermap import supermap
from jcvi.formats.blast import BlastSlow, Blast
from jcvi.formats.sizes import Sizes
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('estsummary', 'provide summary of mapped ESTs with id%, cov% cutoffs'),
        ('rnaseqbench', 'evaluate completeness, contiguity, and chimer of RNAseq'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def estsummary(args):
    """
    %prog estsummary estblastfile estfastafile

    Provide summary on percentage mapped, based on cutoff of id% and cov%. Often
    used in evaluating completeness of genome assembly based on EST mapping.
    `estblastfile` contains the mapping from estfastafile to the genome.
    `estfastafile` contains the FASTA records of the ESTs.
    """
    p = OptionParser(estsummary.__doc__)
    p.add_option("--iden", dest="ident", type="int", default=95,
            help="Cutoff to call mapping [default: %default%]")
    p.add_option("--cov", dest="cov", type="int", default=90,
            help="Cutoff to call mapping [default: %default%]")
    p.add_option("--list", dest="list", default=False, action="store_true",
            help="List the id% and cov% per gene [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    blastfile, fastafile = args
    ident = opts.ident
    cov = opts.cov
    sizes = Sizes(fastafile).mapping
    total = len(sizes)

    covered = 0
    mismatches = 0
    gaps = 0
    alignlen = 0
    queries = set()
    valid_count = 0
    mapped_count = 0
    blast = BlastSlow(blastfile)

    for query, blines in blast.iter_hits():
        blines = list(blines)
        queries.add(query)
        mapped_count += 1

        # per gene report
        this_covered = 0
        this_alignlen = 0
        this_mismatches = 0
        this_gaps = 0

        for b in blines:
            this_covered += abs(b.qstart - b.qstop + 1)
            this_alignlen += b.hitlen
            this_mismatches += b.nmismatch
            this_gaps += b.ngaps

        this_identity = (100 - (this_mismatches + this_gaps) * 100. / this_alignlen)
        this_coverage = this_covered * 100. / sizes[query]
        if opts.list:
            print "{0}\t{1:.1f}\t{2:.1f}".format(query, this_identity, this_coverage)

        if this_identity > ident and this_coverage > cov:
            valid_count += 1

        covered += this_covered
        mismatches += this_mismatches
        gaps += this_gaps
        alignlen += this_alignlen

    print >> sys.stderr, "Cutoff: {0}% ident, {1}% cov".format(ident, cov)
    print >> sys.stderr, "Identity: {0} mismatches, {1} gaps, {2} alignlen".\
            format(mismatches, gaps, alignlen)
    print >> sys.stderr, "Total mapped: {0} ({1:.1f}% of {2})".\
            format(mapped_count, mapped_count * 100. / total, total)
    print >> sys.stderr, "Total valid: {0} ({1:.1f}% of {2})".\
            format(valid_count, valid_count * 100. / total, total)
    print >> sys.stderr, "Overal ident% = {0:.1f}%".\
            format(100 - (mismatches + gaps) * 100. / alignlen)

    queries_combined = sum(sizes[x] for x in queries)
    print >> sys.stderr, "Coverage: {0} covered, {1} total".\
            format(covered, queries_combined)
    print >> sys.stderr, "Coverage = {0:.1f}%".\
            format(covered * 100. / queries_combined)


def rnaseqbench(args):
    """
    %prog rnaseqbench blastfile ref.fasta

    Evaluate de-novo RNA-seq assembly against a reference gene set (same of
    closely related organism). Ideally blatfile needs to be supermap'd.

    Following metric is used (Martin et al. 2010, Rnnotator paper):
    Accuracy: % of contigs share >=95% identity with ref genome (TODO)
    Completeness: % of ref genes covered by contigs to >=80% of their lengths
    Contiguity: % of ref genes covered by a *single* contig >=80% of lengths
    Chimer: % of contigs that contain two or more annotated genes >= 50bp
    """
    p = OptionParser(rnaseqbench.__doc__)
    p.add_option("--cov", dest="cov", type="int", default=80,
            help="Cutoff to call mapping [default: %default%]")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    blastfile, reffasta = args
    cov = opts.cov
    sizes = Sizes(reffasta).mapping
    know_genes = len(sizes)

    querysupermap = blastfile + ".query.supermap"

    # I use extra long names here to prevent cyclic imports
    if not op.exists(querysupermap):
        supermap(blastfile, filter="query")

    blast = Blast(querysupermap)
    ctgs = 0
    chimers = 0
    for ctg, hits in blast.iter_hits():
        bps = defaultdict(int)
        for x in hits:
            bps[x.subject] += abs(x.sstop - x.sstart) + 1

        valid_hits = [x for (x, length) in bps.items() if length >= 100]
        if len(valid_hits) > 1:
            chimers += 1
        ctgs += 1

    print >> sys.stderr, "Chimers: {0} out of {1} ({2:.1f}%).".\
            format(chimers, ctgs, chimers * 100. / ctgs)


if __name__ == '__main__':
    main()
