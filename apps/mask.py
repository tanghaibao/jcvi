#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog input.fasta 

This script pipelines the windowmasker in ncbi-blast+.
The masked fasta file will have an appended suffix of .mask with all the
low-complexity nucleotides soft-masked (to lower case)
"""
import os
import sys

from jcvi.apps.base import sh, is_current_file, debug
debug()


def wm_mk_counts(genome_file):
    ext = ".counts"
    cmd = "windowmasker -in %(genome_file)s -mk_counts " \
          "-out %(genome_file)s%(ext)s" % locals() 

    if not is_current_file(genome_file + ext, genome_file):
        sh(cmd)

def wm_mk_masks(genome_file):
    ext = ".mask"
    outfmt = "fasta"
    cmd = "windowmasker -in %(genome_file)s -ustat %(genome_file)s.counts " \
          "-outfmt %(outfmt)s -dust T " \
          "-out %(genome_file)s%(ext)s" % locals()

    if not is_current_file(genome_file + ext, genome_file):
        sh(cmd)


def main():

    from optparse import OptionParser
    
    p = OptionParser(__doc__)
    opts, args = p.parse_args()
    
    argc = len(args)
    if argc != 1:
        sys.exit(p.print_help())

    genome_file = args[0] 

    # entire pipeline
    wm_mk_counts(genome_file)
    wm_mk_masks(genome_file)


if __name__ == '__main__':
    main()
