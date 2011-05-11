#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Prepare input files for Celera Assembler, dispatch based on file suffix::

*.fasta: convert-fasta-to-v2.pl
*.sff: sffToCA
*.fastq: fastqToCA
"""

import os.path as op
import sys
import logging

from optparse import OptionParser

from Bio import SeqIO

from jcvi.formats.fasta import get_qual, iter_fasta_qual, write_fasta_qual
from jcvi.apps.softlink import get_abs_path
from jcvi.apps.base import ActionDispatcher, sh, set_grid, debug
debug()

CAPATH = "~/bin/Linux-amd64/bin/"


def main():

    actions = (
        ('fasta', 'convert fasta to frg file'),
        ('sff', 'convert 454 reads to frg file'),
        ('fastq', 'convert Illumina reads to frg file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def make_qual(fastafile, defaultqual=20):
    """
    Make a qualfile with default qual value if not available
    """
    assert op.exists(fastafile)

    qualfile = get_qual(fastafile)
    if qualfile is None:
        qualfile = get_qual(fastafile, check=False)
        qualhandle = open(qualfile, "w")

        for rec in iter_fasta_qual(fastafile, None, defaultqual=defaultqual):
            write_fasta_qual(rec, None, qualhandle)

        logging.debug("write qual values to file `{0}`".format(qualfile))
        qualhandle.close()

    return qualfile


def make_matepairs(fastafile):
    """
    Assumes the mates are adjacent sequence records
    """
    assert op.exists(fastafile)

    matefile = fastafile.rsplit(".", 1)[0] + ".mates"
    if op.exists(matefile):
        logging.debug("matepairs file `{0}` found".format(matefile))
    else:
        logging.debug("parsing matepairs from `{0}`".format(fastafile))
        matefw = open(matefile, "w")
        it = SeqIO.parse(fastafile, "fasta")
        for fwd, rev in zip(it, it):
            print >> matefw, "{0}\t{1}".format(fwd.id, rev.id)

        matefw.close()

    return matefile


def add_size_option(p):
    p.add_option("-s", dest="size", default=0, type="int",
            help="insert has mean size of [default: %default] " + \
                 "stddev is assumed to be 25% around mean size")


get_mean_sv = lambda size: (size, size / 4)


def fasta(args):
    """
    %prog fasta fastafile

    Convert reads formatted as FASTA file, and convert to CA frg file. If .qual
    file is found, then use it, otherwise just make a fake qual file. Mates are
    assumed as adjacent sequence records (i.e. /1, /2, /1, /2 ...) unless a
    matefile is given.
    """
    p = OptionParser(fasta.__doc__)
    p.add_option("-m", dest="matefile", default=None,
            help="matepairs file")
    set_grid(p)
    add_size_option(p)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    grid = opts.grid

    fastafile, = args
    libname = op.basename(fastafile).split(".")[0]
    frgfile = fastafile.rsplit(".", 1)[0] + ".frg"

    mated = (opts.size != 0)
    mean, sv = get_mean_sv(opts.size)

    qualfile = make_qual(fastafile)
    if mated:
        if opts.matefile:
            matefile = opts.matefile
            assert op.exists(matefile)
        else:
            matefile = make_matepairs(fastafile)

    cmd = CAPATH + "convert-fasta-to-v2.pl -l {0} -s {1} -q {2} ".\
            format(libname, fastafile, qualfile)
    if mated:
        cmd += "-mean {0} -stddev {1} -m {2} ".format(mean, sv, matefile)

    sh(cmd, grid=grid, outfile=frgfile)


def sff(args):
    """
    %prog sff sffiles

    Convert reads formatted as 454 SFF file, and convert to CA frg file.
    """
    p = OptionParser(sff.__doc__)
    p.add_option("-o", dest="output", default="out",
            help="output frg filename")
    set_grid(p)
    add_size_option(p)

    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(p.print_help())

    grid = opts.grid

    sffiles = args
    libname = opts.output

    mated = (opts.size != 0)
    mean, sv = get_mean_sv(opts.size)

    cmd = CAPATH + "sffToCA -libraryname {0} -output {0} ".format(libname)
    cmd += "-clear 454 -trim chop "
    if mated:
        cmd += "-linker titanium -insertsize {0} {1} ".format(mean, sv)

    cmd += " ".join(sffiles)

    sh(cmd, grid=grid)


def fastq(args):
    """
    %prog fastq fastqfile

    Convert reads formatted as FASTQ file, and convert to CA frg file.
    """
    p = OptionParser(fastq.__doc__)
    p.add_option("--sanger", dest="sanger", default=False, action="store_true",
            help="Are the qv sanger encodings? [default: %default]")
    p.add_option("--outtie", dest="outtie", default=False, action="store_true",
            help="Are these outie reads? [default: %default]")
    add_size_option(p)

    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(p.print_help())

    fastqfile = get_abs_path(args[0])
    if len(args) == 2:
        fastqfile2 = get_abs_path(args[1])

    libname = op.basename(fastqfile).split(".")[0]
    frgfile = libname + ".frg"

    mated = (opts.size != 0)
    mean, sv = get_mean_sv(opts.size)

    cmd = CAPATH + "fastqToCA -libraryname {0} -fastq {1}".\
            format(libname, fastqfile)
    if mated:
        assert len(args) == 2, "you need two fastq file for mated library"
        cmd += ",{0} -insertsize {1} {2} ".format(fastqfile2, mean, sv)
    if opts.sanger:
        cmd += " -type sanger "
    if opts.outtie:
        cmd += " -outtie "

    sh(cmd, outfile=frgfile)


if __name__ == '__main__':
    main()
