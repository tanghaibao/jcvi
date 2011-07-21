#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run through NCBI vecscreen on a local machine.
"""

import sys
import shutil
import logging

from optparse import OptionParser

from jcvi.utils.cbook import depends 
from jcvi.utils.range import range_merge
from jcvi.formats.fasta import tidy
from jcvi.formats.blast import BlastLine
from jcvi.formats.base import must_open
from jcvi.apps.base import ActionDispatcher, debug, download, sh
debug()


def main():

    actions = (
        ('blast', 'run BLASTN against UniVec'),
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
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    outfastafile = fastafile.rsplit(".", 1)[0] + ".masked.fasta"
    vecbedfile = blast([fastafile])
    ecolifile = "Ecoli.fasta"
    assert op.exists(ecolifile)
    ecolibedfile = blast([fastafile, "--db={0}".format(ecolifile), "--blat"])

    cmd = "cat {0} {1}".format(vecbedfile, ecolibedfile)
    cmd += " | mergeBed -nms -d 100 -i stdin"
    cmd += " | maskFastaFromBed -fi {0} -bed stdin -fo {1}".\
            format(fastafile, outfastafile)
    sh(cmd)

    tidy([outfastafile])


@depends
def run_formatdb(infile=None, outfile=None):
    cmd = "formatdb -i {0} -p F".format(infile)
    sh(cmd)


@depends
def run_blastall(infile=None, outfile=None, db="UniVec_Core"):
    """
    BLASTN parameters reference:
    http://www.ncbi.nlm.nih.gov/VecScreen/VecScreen_docs.html
    """
    cmd = 'blastall -p blastn -i {0}'.format(infile) 
    cmd += ' -d {0} -q -5 -G 3 -E 3 -F "m D"'.format(db)
    cmd += ' -e 0.01 -Y 1.75e12 -m 8 -o {0} -a 8'.format(outfile)
    sh(cmd)


@depends
def run_blat(infile=None, outfile=None, db="UniVec_Core"):
    cmd = 'blat {0} {1} -out=blast8 {2}'.format(db, infile, outfile)
    sh(cmd)

    blatfile = outfile
    filtered_blatfile = outfile + ".P95L50"
    run_blast_filter(infile=blatfile, outfile=filtered_blatfile)
    shutil.move(filtered_blatfile, blatfile)


def run_blast_filter(infile=None, outfile=None):
    from jcvi.formats.blast import filter

    logging.debug("Filter BLAST result (pctid=95, hitlen=50)")
    filter([infile, "--pctid=95", "--hitlen=50"])


def blast(args):
    """
    %prog blast fastafile 

    Run BLASTN against database (default is UniVec_Core).  Output .bed format
    on the vector/contaminant ranges.
    """
    p = OptionParser(blast.__doc__)
    p.add_option("--dist", dest="dist", default=100, type="int",
            help="Merge adjacent HSPs separated by [default: %default]")
    p.add_option("--db", dest="db", default=None,
            help="Use a different database rather than UniVec_Core")
    p.add_option("--blat", dest="blat", default=False, action="store_true",
            help="Use BLAT instead of BLAST [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    fastaprefix = fastafile.split(".", 1)[0]

    univec = opts.db or download("ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core")
    univecdb = univec + ".nin"
    uniprefix = univec.split(".", 1)[0]
    run_formatdb(infile=univec, outfile=univecdb)

    fastablast = fastaprefix + ".{0}.blast".format(uniprefix)

    prog = run_blat if opts.blat else run_blastall
    prog(infile=fastafile, outfile=fastablast, db=univec)

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
