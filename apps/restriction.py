#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Procedure to cut genome using restriction enzymes.
"""

import sys
import logging

from optparse import OptionParser

from Bio.Restriction.Restriction import AllEnzymes, Analysis

from jcvi.formats.fasta import Fasta, SeqRecord, SeqIO
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('fragment', 'extract upstream and downstream seq of particular RE'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fragment(args):
    """
    %prog fragment fastafile enzyme

    Cut the fastafile using the specified enzyme, and grab upstream and
    downstream nucleotide sequence along with the cut site.
    """
    p = OptionParser(fragment.__doc__)
    p.add_option("--flank", default=10, type="int",
            help="Extract flanking bases of the cut sites [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, enzyme = args
    flank = opts.flank
    assert flank > 0

    assert enzyme in set(str(x) for x in AllEnzymes)
    fragfastafile = fastafile.split(".")[0] + \
        ".{0}.flank{1}.fasta".format(enzyme, flank)
    enzyme = [x for x in AllEnzymes if str(x) == enzyme][0]
    fw = open(fragfastafile, "w")

    f = Fasta(fastafile, lazy=True)
    for name, rec in f.iteritems_ordered():
        a = Analysis([enzyme], rec.seq)
        sites = a.full()[enzyme]
        size = len(rec)
        for s in sites:
            newid = "{0}:{1}".format(name, s)
            left = max(s - flank, 0)
            right = min(s + flank, size)
            frag = rec.seq[left:right]
            frag = frag.strip("Nn")
            newrec = SeqRecord(frag, id=newid, description="")
            SeqIO.write([newrec], fw, "fasta")

    logging.debug("Fragments written to `{0}`.".format(fragfastafile))


if __name__ == '__main__':
    main()
