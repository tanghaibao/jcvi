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

ECOLI_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/019/425/GCF_000019425.1_ASM1942v1/GCF_000019425.1_ASM1942v1_genomic.fna.gz"
UNIVEC_URL = "ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core"


def main():

    actions = (("mask", "mask the contaminants"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def is_internet_file(url):
    """Return if url starts with http://, https://, or ftp://.

    Args:
        url (str): URL of the link
    """
    return (
        url.startswith("http://")
        or url.startswith("https://")
        or url.startswith("ftp://")
    )


def mask(args):
    """
    %prog mask fastafile

    Mask the contaminants. By default, this will compare against UniVec_Core and
    Ecoli.fasta. Merge the contaminant results, and use `maskFastaFromBed`. Can
    perform FASTA tidy if requested.
    """
    p = OptionParser(mask.__doc__)
    p.add_option(
        "--db",
        default=ECOLI_URL,
        help="Contaminant db other than Ecoli K12, will download if file starts with http://, https://, or ftp://",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    db = opts.db
    assert op.exists(fastafile)

    outfastafile = fastafile.rsplit(".", 1)[0] + ".masked.fasta"
    vecbedfile = blast([fastafile])
    ecolifile = (
        download(db, filename="Ecoli.fasta", handle_gzip=True)
        if is_internet_file(db)
        else db
    )
    assert op.exists(ecolifile)
    ecolibedfile = blast([fastafile, "--db={0}".format(ecolifile)])

    cmd = "cat {0} {1}".format(vecbedfile, ecolibedfile)
    cmd += " | sort -k1,1 -k2,2n"
    cmd += " | mergeBed -c 4 -o distinct -d 100 -i stdin"
    cmd += " | maskFastaFromBed -fi {0} -bed stdin -fo {1}".format(
        fastafile, outfastafile
    )
    sh(cmd)

    return tidy([outfastafile])


def blast(args):
    """
    %prog blast fastafile

    Run BLASTN against database (default is UniVec_Core).  Output .bed format
    on the vector/contaminant ranges.
    """
    p = OptionParser(blast.__doc__)
    p.add_option(
        "--dist",
        default=100,
        type="int",
        help="Merge adjacent HSPs separated by",
    )
    p.add_option("--db", help="Use a different database rather than UniVec_Core")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    fastaprefix = fastafile.split(".", 1)[0]

    univec = opts.db or download(UNIVEC_URL)
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
        print("\t".join(str(x) for x in (seqid, start - 1, end, uniprefix)), file=fw)

    return bedfile


if __name__ == "__main__":
    main()
