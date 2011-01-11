#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog input.fasta 

This script pipelines the windowmasker in ncbi-blast+.
The masked fasta file will have an appended suffix of .mask with all the
low-complexity nucleotides soft-masked (to lower case)
"""

import os
import os.path as op
import sys

from jcvi.apps.base import sh, is_current_file, debug
debug()


def wm_mk_counts(genomefile):
    outfile = "%s.counts" % genomefile
    cmd = "windowmasker -in %(genomefile)s -mk_counts " \
          "-out %(outfile)s" % locals() 

    if not is_current_file(outfile, genomefile):
        sh(cmd)

def wm_mk_masks(genomefile):
    outfile = "%s.masked%s" % op.splitext(genomefile)
    cmd = "windowmasker -in %(genomefile)s -ustat %(genomefile)s.counts " \
          "-outfmt fasta -dust T " \
          "-out %(outfile)s" % locals()

    if not is_current_file(outfile, genomefile):
        sh(cmd)


def main():

    from optparse import OptionParser
    
    p = OptionParser(__doc__)
    opts, args = p.parse_args()
    
    argc = len(args)
    if argc != 1:
        sys.exit(p.print_help())

    genomefile = args[0] 

    # entire pipeline
    wm_mk_counts(genomefile)
    wm_mk_masks(genomefile)


if __name__ == '__main__':
    main()
