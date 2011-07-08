#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
"""

import os
import os.path as op
import sys

from optparse import OptionParser

from jcvi.formats.fasta import Fasta
from jcvi.utils.cbook import percentage
from jcvi.apps.base import ActionDispatcher, debug, set_grid, sh, is_newer_file
debug()


def wm_mk_counts(genomefile):
    outfile = "%s.counts" % genomefile
    cmd = "windowmasker -in %(genomefile)s -mk_counts " \
          "-out %(outfile)s" % locals()

    if not is_newer_file(outfile, genomefile):
        sh(cmd)

    return outfile


def wm_mk_masks(genomefile):
    outfile = "%s.masked%s" % op.splitext(genomefile)
    cmd = "windowmasker -in %(genomefile)s -ustat %(genomefile)s.counts " \
          "-outfmt fasta -dust T " \
          "-out %(outfile)s" % locals()

    sh(cmd)

    return outfile


def hardmask(fastafile):
    cmd = "maskOutFa {0} hard {0}".format(fastafile)
    sh(cmd)


def main():

    actions = (
        ('mask', 'use windowmasker to mask low-complexity bases'),
        ('summary', 'report the number of bases and sequences masked'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def summary(args):
    """
    %prog summary fastafile

    Report the number of bases and sequences masked.
    """
    p = OptionParser(summary.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    f = Fasta(fastafile, index=False)

    halfmaskedseqs = set()
    allmasked = 0
    allbases = 0
    cutoff = 50
    others = 0
    for key, seq in f.iteritems():
        masked = 0
        for base in seq:
            if base not in "AGCT":
                masked += 1
        seqlen = len(seq)
        if masked * 100. / seqlen > cutoff:
            halfmaskedseqs.add(key)
        allmasked += masked
        allbases += seqlen

    seqnum = len(f)
    maskedseqnum = len(halfmaskedseqs)

    print >> sys.stderr, "Total masked bases: {0}".\
            format(percentage(allmasked, allbases))
    print >> sys.stderr, "Total masked sequences (contain > {0}% masked): {1}".\
            format(cutoff, percentage(maskedseqnum, seqnum))


def mask(args):
    """
    %prog mask fastafile

    This script pipelines the windowmasker in NCBI BLAST+.  Masked fasta file
    will have an appended suffix of .mask with all the low-complexity bases masked
    (default to lower case, set --hard for hardmasking).
    """
    p = OptionParser(mask.__doc__)
    p.add_option("--hard", dest="hard", default=False, action="store_true",
            help="Hard mask the low-complexity bases [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    genomefile, = args

    # entire pipeline
    wm_mk_counts(genomefile)
    maskedfastafile = wm_mk_masks(genomefile)

    if opts.hard:
        hardmask(maskedfastafile)


if __name__ == '__main__':
    main()
