"""
%prog blast_file cds_file bed_file [options]

Find tandem gene clusters that are separated by N genes, based on filtered
blast_file by enforcing the alignments between any two genes are at least 50% of
either genes.
"""

import sys
import logging

from optparse import OptionParser

from jcvi.formats.blast import BlastLine
from jcvi.formats.fasta import Fasta
from jcvi.formats.bed import Bed
from jcvi.utils.grouper import Grouper
from jcvi.utils.cbook import gene_name
from jcvi.apps.base import debug
debug()


def main(blast_file, cds_file, bed_file, N=3):

    # get the sizes for the CDS first
    f = Fasta(cds_file)
    sizes = dict(f.itersizes())
    
    # retrieve the locations
    bed = Bed(bed_file).order

    # filter the blast file
    g = Grouper()
    fp = open(blast_file)
    for row in fp:
        b = BlastLine(row)
        query_len = sizes[b.query]
        subject_len = sizes[b.subject]
        if b.hitlen < min(query_len, subject_len) / 2:
            continue
        
        query, subject = gene_name(b.query), gene_name(b.subject)
        qi, q = bed[query]
        si, s = bed[subject]

        if q.seqid==s.seqid and abs(qi-si) <= N:
            g.join(query, subject)

    # dump the grouper
    ngenes, nfamilies = 0, 0
    families = []
    for group in sorted(g):
        if len(group) >= 2:
            print ",".join(sorted(group))
            ngenes += len(group)
            nfamilies += 1 
            families.append(sorted(group))

    longest_family = max(families, key=lambda x: len(x))
        
    # generate reports
    print >>sys.stderr, "Proximal paralogues (dist=%d):" % N
    print >>sys.stderr, "Total of %d genes in %d families" % (ngenes, nfamilies)
    print >>sys.stderr, "Longest families (%d): %s" % (len(longest_family),
        ",".join(longest_family))


if __name__ == '__main__':

    p = OptionParser(__doc__)
    p.add_option("--tandem_Nmax", dest="tandem_Nmax", type="int", default=3,
               help="merge tandem genes within distance [default: %default]")
    
    (opts, args) = p.parse_args()
    
    if len(args) != 3:
        sys.exit(p.print_help())

    blast_file, cds_file, bed_file = args
    N = opts.tandem_Nmax

    main(blast_file, cds_file, bed_file, N=N)

