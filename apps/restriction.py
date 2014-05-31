#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Procedure to cut genome using restriction enzymes.
"""

import sys
import logging

from Bio.Restriction.Restriction import AllEnzymes, Analysis

from jcvi.formats.fasta import Fasta, SeqRecord, SeqIO
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('fragment', 'extract upstream and downstream seq of particular RE'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def extract_full(rec, sites, flank, fw):
    """
    Full extraction of seq flanking the sites.
    """
    for s in sites:
        newid = "{0}:{1}".format(rec.name, s)
        left = max(s - flank, 0)
        right = min(s + flank, len(rec))
        frag = rec.seq[left:right].strip("Nn")
        newrec = SeqRecord(frag, id=newid, description="")
        SeqIO.write([newrec], fw, "fasta")


def extract_ends(rec, sites, flank, fw, maxfragsize=800):
    """
    Extraction of ends of fragments above certain size.
    """
    nsites = len(sites)
    size = len(rec)
    for i, s in enumerate(sites):
        newid = "{0}:{1}".format(rec.name, s)
        recs = []

        if i == 0 or s - sites[i - 1] <= maxfragsize:
            newidL = newid + "L"
            left = max(s - flank, 0)
            right = s
            frag = rec.seq[left:right].strip("Nn")
            recL = SeqRecord(frag, id=newidL, description="")
            if i == 0 and s > maxfragsize:  # Contig L-end
                pass
            else:
                recs.append(recL)

        if i == nsites - 1 or sites[i + 1] - s <= maxfragsize:
            newidR = newid + "R"
            left = s
            right = min(s + flank, size)
            frag = rec.seq[left:right].strip("Nn")
            recR = SeqRecord(frag, id=newidR, description="")
            if i == nsites - 1 and size - s > maxfragsize:  # Contig R-end
                pass
            else:
                recs.append(recR)

        SeqIO.write(recs, fw, "fasta")


def fragment(args):
    """
    %prog fragment fastafile enzyme

    Cut the fastafile using the specified enzyme, and grab upstream and
    downstream nucleotide sequence along with the cut site. In this case, the
    sequences extracted are:

                |- PstI
    ============|===========
            (-------)

    Sometimes we need to limit the size of the restriction fragments, for
    example the GBS protocol does not allow fragments larger than 800bp.

           |-PstI        |- PstI              |- PstI
    ~~~====|=============|==========~~~~~~~===|============
           (---)     (---)

    In this case, the second fragment is longer than 800bp, therefore the two
    ends are NOT extracted, as in the first fragment.
    """
    p = OptionParser(fragment.__doc__)
    p.add_option("--flank", default=150, type="int",
            help="Extract flanking bases of the cut sites [default: %default]")
    p.add_option("--full", default=False, action="store_true",
            help="The full extraction mode [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, enzyme = args
    flank = opts.flank
    assert flank > 0
    extract = extract_full if opts.full else extract_ends
    tag = "full" if opts.full else "ends"

    assert enzyme in set(str(x) for x in AllEnzymes)
    fragfastafile = fastafile.split(".")[0] + \
        ".{0}.flank{1}.{2}.fasta".format(enzyme, flank, tag)
    enzyme = [x for x in AllEnzymes if str(x) == enzyme][0]

    f = Fasta(fastafile, lazy=True)
    fw = open(fragfastafile, "w")
    for name, rec in f.iteritems_ordered():
        a = Analysis([enzyme], rec.seq)
        sites = a.full()[enzyme]
        extract(rec, sites, flank, fw)

    logging.debug("Fragments written to `{0}`.".format(fragfastafile))


if __name__ == '__main__':
    main()
