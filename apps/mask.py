#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog.py input.fasta <species_name>

This script pipelines the windowmasker in ncbi-blast+.
Species_name is optional, when added, a blast database
with descriptions including the species_name
"""
import os
from subprocess import Popen
import sys

import logging
logging.basicConfig(level=logging.DEBUG)

from jcvi.apps.base import sh, is_current_file


def wm_mk_counts(genome_file):
    ext = ".counts"
    cmd = "windowmasker -in %(genome_file)s -mk_counts -dust T " \
          "-out %(genome_file)s%(ext)s" % locals()

    if not is_current_file(genome_file + ext, genome_file):
        sh(cmd)

def wm_mk_masks(genome_file, mask_fasta=True):
    ext = ".mask" if mask_fasta else ".asnb"
    outfmt = "fasta" if mask_fasta else "maskinfo_asn1_bin"
    cmd = "windowmasker -in %(genome_file)s -ustat %(genome_file)s.counts " \
          "-outfmt %(outfmt)s " \
          "-out %(genome_file)s%(ext)s" % locals()

    if not is_current_file(genome_file + ext, genome_file):
        sh(cmd)

def wm_mk_database(genome_file, species):
    ext = ".nin"
    cmd = "makeblastdb -in %(genome_file)s -dbtype nucl " \
          "-mask_data %(genome_file)s.asnb -out %(genome_file)s " \
          " -title \"%(species)s genome, wm masked \"" % locals()

    if not is_current_file(genome_file + ext, genome_file):
        sh(cmd)


def main():

    from optparse import OptionParser
    
    p = OptionParser(__doc__)
    opts, args = p.parse_args()
    
    argc = len(args)
    if argc not in (1, 2):
        sys.exit(p.print_help())

    genome_file = args[0] 
    if argc==2:
        species = args[1]

    # entire pipeline
    wm_mk_counts(genome_file)
    wm_mk_masks(genome_file)

    if argc==2:
        wm_mk_database(genome_file, species)


if __name__ == '__main__':
    main()
