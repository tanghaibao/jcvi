#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run through NCBI vecscreen on a local machine.
"""

import os.path as op
import sys

from jcvi.utils.range import range_merge
from jcvi.formats.fasta import tidy
from jcvi.formats.blast import BlastLine
from jcvi.formats.base import must_open
from jcvi.apps.align import run_vecscreen, run_megablast
from jcvi.apps.base import OptionParser, ActionDispatcher, download, sh


def main():

    actions = (
        ('mask', 'mask the contaminants'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def mask(args):
    """
    %prog mask fastafile

    Mask the contaminants. By default, this will compare against UniVec_Core and
    Ecoli.fasta. Merge the contaminant results, and use `maskFastaFromBed`. Can
    perform FASTA tidy if requested.
    """
    p = OptionParser(mask.__doc__)
    p.add_option("--db",
                 help="Contaminant db other than Ecoli K12 [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    assert op.exists(fastafile)

    outfastafile = fastafile.rsplit(".", 1)[0] + ".masked.fasta"
    vecbedfile = blast([fastafile])
    ecoliurl = \
    "ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid58979/NC_010473.fna"
    ecolifile = opts.db or download(ecoliurl, filename="Ecoli.fasta")
    assert op.exists(ecolifile)
    ecolibedfile = blast([fastafile, "--db={0}".format(ecolifile)])

    cmd = "cat {0} {1}".format(vecbedfile, ecolibedfile)
    cmd += " | mergeBed -nms -d 100 -i stdin"
    cmd += " | maskFastaFromBed -fi {0} -bed stdin -fo {1}".\
            format(fastafile, outfastafile)
    sh(cmd)

    return tidy([outfastafile])


def blast(args):
    """
    %prog blast fastafile

    Run BLASTN against database (default is UniVec_Core).  Output .bed format
    on the vector/contaminant ranges.
    """
    p = OptionParser(blast.__doc__)
    p.add_option("--dist", default=100, type="int",
            help="Merge adjacent HSPs separated by [default: %default]")
    p.add_option("--db",
            help="Use a different database rather than UniVec_Core")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    fastaprefix = fastafile.split(".", 1)[0]

    univec = opts.db or download("ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core")
    uniprefix = univec.split(".", 1)[0]

    fastablast = fastaprefix + ".{0}.blast".format(uniprefix)

    prog = run_megablast if opts.db else run_vecscreen
    prog(infile=fastafile, outfile=fastablast, db=univec, pctid=95, hitlen=50)

    fp = open(fastablast)
    ranges = []
    for row in fp:
        b = BlastLine(row)
        ranges.append((b.query, b.qstart, b.qstop))

    merged_ranges = range_merge(ranges, dist=opts.dist)
    bedfile = fastaprefix + ".{0}.bed".format(uniprefix)
    fw = must_open(bedfile, "w")
    for seqid, start, end in merged_ranges:
        print >> fw, "\t".join(str(x) for x in (seqid, start - 1, end, uniprefix))

    return bedfile


if __name__ == '__main__':
    main()
