#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Mask low complexity regions in the genome.
"""
from __future__ import print_function

import os.path as op
import sys

from jcvi.formats.fasta import Fasta
from jcvi.utils.cbook import depends, percentage
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


@depends
def wm_mk_counts(infile=None, outfile=None):
    cmd = "windowmasker -in {0} -mk_counts".format(infile)
    cmd += " -out {0}".format(outfile)
    sh(cmd)


@depends
def wm_mk_masks(infile=None, outfile=None, genomefile=None):
    cmd = "windowmasker -in {0} -ustat {1}".format(genomefile, infile)
    cmd +=  " -outfmt fasta -dust T -out {0}".format(outfile)
    sh(cmd)


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

    print("Total masked bases: {0}".\
            format(percentage(allmasked, allbases)), file=sys.stderr)
    print("Total masked sequences (contain > {0}% masked): {1}".\
            format(cutoff, percentage(maskedseqnum, seqnum)), file=sys.stderr)


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
    countsfile = genomefile + ".counts"
    wm_mk_counts(infile=genomefile, outfile=countsfile)

    maskedfastafile = "%s.masked%s" % op.splitext(genomefile)
    wm_mk_masks(infile=countsfile, outfile=maskedfastafile, genomefile=genomefile)

    if opts.hard:
        hardmask(maskedfastafile)


if __name__ == '__main__':
    main()
