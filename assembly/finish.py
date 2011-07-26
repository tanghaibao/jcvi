#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Finishing pipeline, starting with a phase1/2 BAC. The pipeline ideally should
include the following components

+ BLAST against the Illumina contigs to fish out additional seqs
+ Use minimus2 to combine the contigs through overlaps
+ Map the mates to the contigs and perform scaffolding
+ Base corrections using Illumina reads
"""

import os
import os.path as op
import sys
import shutil

from optparse import OptionParser

from jcvi.formats.fasta import gaps
from jcvi.utils.cbook import depends
from jcvi.apps.base import run_formatdb, run_blast_filter
from jcvi.apps.base import ActionDispatcher, debug, sh
debug()


def main():

    actions = (
        ('overlap', 'Build larger contig set by fishing additional seqs'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


@depends
def run_gapsplit(infile=None, outfile=None):
    gaps([infile, "--split"])
    return outfile


@depends
def run_megablast(infile=None, outfile=None, db=None):
    cmd = "megablast -i {0} -d {1}".format(infile, db)
    cmd += " -e 0.001 -m 8 -o {0}".format(outfile)
    sh(cmd)

    blastfile = outfile
    filtered_blastfile = outfile + ".P98L100"
    run_blast_filter(infile=blastfile, outfile=filtered_blastfile,
            pctid=98, hitlen=100)
    shutil.move(filtered_blastfile, blastfile)


def overlap(args):
    """
    %prog overlap ctgfasta poolfasta

    Fish out the sequences in `poolfasta` that overlap with `ctgfasta`.
    Mix and combine using `minimus2`.
    """
    p = OptionParser(overlap.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ctgfasta, poolfasta = args
    prefix = ctgfasta.split(".")[0]
    splitctgfasta = ctgfasta.rsplit(".", 1)[0] + ".split.fasta"
    ctgfasta = run_gapsplit(infile=ctgfasta, outfile=splitctgfasta)

    poolfastadb = poolfasta + ".nin"
    run_formatdb(infile=poolfasta, outfile=poolfastadb)

    blastfile = ctgfasta + ".blast"
    run_megablast(infile=ctgfasta, outfile=blastfile, db=poolfasta)

    cmd = "cut -f2 {0} | sort -u > {1}.ids".format(blastfile, prefix)
    sh(cmd)

    cmd = "faSomeRecords {0} {1}.ids {1}.ids.fasta".format(poolfasta, prefix)
    sh(cmd)

    cmd = "cat {0} {1}.ids.fasta > {1}.merged.fasta".format(ctgfasta, prefix)
    sh(cmd)

    cmd = "toAmos -s {0}.merged.fasta -o {0}.afg".format(prefix)
    sh(cmd)

    cmd = "minimus2 {0} -D REFCOUNT=0".format(prefix)
    sh(cmd)


if __name__ == '__main__':
    main()
