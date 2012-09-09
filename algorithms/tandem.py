"""
%prog blast_file cds_file bed_file [options]

Find tandem gene clusters that are separated by N genes, based on filtered
blast_file by enforcing alignments between any two genes at least 50%
(or user specified value) of either gene.

pep_file can also be used in same manner.
"""

import sys
import logging

from optparse import OptionParser

from jcvi.formats.blast import BlastLine
from jcvi.formats.fasta import Fasta
from jcvi.formats.bed import Bed
from jcvi.formats.base import must_open
from jcvi.utils.grouper import Grouper
from jcvi.utils.cbook import gene_name
from jcvi.apps.base import debug, set_outfile
debug()


def main(blast_file, cds_file, bed_file, N=3, P=50, is_self=True, \
    strip_name=".", ofile=sys.stderr):

    # get the sizes for the CDS first
    f = Fasta(cds_file)
    sizes = dict(f.itersizes())

    # retrieve the locations
    bed = Bed(bed_file)
    order = bed.order

    if is_self:
        # filter the blast file
        g = Grouper()
        fp = open(blast_file)
        for row in fp:
            b = BlastLine(row)
            query_len = sizes[b.query]
            subject_len = sizes[b.subject]
            if b.hitlen < min(query_len, subject_len)*P/100.:
                continue

            query = gene_name(b.query, strip_name)
            subject = gene_name(b.subject, strip_name)
            qi, q = order[query]
            si, s = order[subject]

            if q.seqid == s.seqid and abs(qi - si) <= N:
                g.join(query, subject)

    else:
        homologs = Grouper()
        fp = open(blast_file)
        for row in fp:
            b = BlastLine(row)
            query_len = sizes[b.query]
            subject_len = sizes[b.subject]
            if b.hitlen < min(query_len, subject_len)*P/100.:
                continue

            query = gene_name(b.query, strip_name)
            subject = gene_name(b.subject, strip_name)
            homologs.join(query, subject)

        g = Grouper()
        for i, atom in enumerate(bed):
            for x in range(1, N+1):
                if all([i-x >= 0, bed[i-x].seqid == atom.seqid, \
                    homologs.joined(bed[i-x].accn, atom.accn)]):
                    g.join(bed[i-x].accn, atom.accn)

    # dump the grouper
    fw = must_open(ofile, "w")
    ngenes, nfamilies = 0, 0
    families = []
    for group in sorted(g):
        if len(group) >= 2:
            print >>fw, ",".join(sorted(group))
            ngenes += len(group)
            nfamilies += 1
            families.append(sorted(group))

    longest_family = max(families, key=lambda x: len(x))

    # generate reports
    print >>sys.stderr, "Proximal paralogues (dist=%d):" % N
    print >>sys.stderr, "Total %d genes in %d families" % (ngenes, nfamilies)
    print >>sys.stderr, "Longest families (%d): %s" % (len(longest_family),
        ",".join(longest_family))

    return families


if __name__ == '__main__':

    p = OptionParser(__doc__)
    p.add_option("--tandem_Nmax", dest="tandem_Nmax", type="int", default=3,
               help="merge tandem genes within distance [default: %default]")
    p.add_option("--percent_overlap", type="int", default=50,
               help="tandem genes have >=x% aligned sequence, x=0-100 \
               [default: %default]")
    p.add_option("--not_self", default=False, action="store_true",
                 help="provided is not self blast file [default: %default]")
    p.add_option("--strip_gene_name", dest="sep", type="string", default=".",
               help="strip alternative splicing. Use None for no stripping. \
               [default: %default]")
    set_outfile(p)

    (opts, args) = p.parse_args()

    if len(args) != 3:
        sys.exit(p.print_help())

    blast_file, cds_file, bed_file = args
    N = opts.tandem_Nmax
    P = opts.percent_overlap
    is_self = not opts.not_self
    sep = opts.sep
    ofile = opts.outfile

    main(blast_file, cds_file, bed_file, N=N, P=P, is_self=is_self, \
        strip_name=sep, ofile=ofile)
