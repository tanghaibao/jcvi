#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Processing fastq files
"""

import os.path as op
import sys

from optparse import OptionParser

from Bio import SeqIO
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('pair', 'pair up two fastq files and combine pairs'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def pair(args):
    """
    %prog pair 1.fastq 2.fastq frags.fastq pairs.fastq 

    pair up the records in 1.fastq and 2.fastq, pairs are indicated by trailing
    "/1" and "/2". If using raw sequences, this is trivial, since we can just
    iterate one at a time for both files; however if two files are not matching,
    (e.g. due to trimming), we need a third fastq that provides the order.
    """
    p = OptionParser(pair.__doc__)
    p.add_option("-Q", dest="infastq", default="illumina", 
            help="input fastq [default: %default]")
    p.add_option("-q", dest="outfastq", default="sanger",
            help="output fastq format [default: %default]")
    p.add_option("-r", dest="ref", default=None,
            help="a reference fastq that provides order")
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(p.print_help())

    afastq, bfastq, frags, pairs = args

    assert op.exists(afastq) and op.exists(bfastq)
    assert not op.exist(frags) and not op.exist(pairs)

    fastqfmt = lambda x: "fastq" if x=="sanger" else "fastq-illumina"
    infmt = fastqfmt(opts.infastq)
    outfmt = fastqfmt(opts.outfastq)
    ref = opts.ref

    ah = SeqIO.parse(afastq, infmt)
    bh = SeqIO.parse(bfastq, infmt)
    if ref: rh = SeqIO.parse(ref, infmt)

    strip_name = lambda x: x.rsplit("/", 1)[0]

    if ref:
        for rec in rh:
            name = strip_name(rec.id)
    else:
        # TODO: unimplemented
        pass


if __name__ == '__main__':
    main()
