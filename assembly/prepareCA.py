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
from jcvi.apps.base import ActionDispatcher, sh, debug
debug()


def main():

    actions = (
        ('fasta', 'convert fasta to frg file'),
        ('sff', 'convert 454 reads to frg file'),
        ('fastq', 'convert Illumina reads to frg file'),
        ('trim', 'CA6.1 cannot handle long fastq seq and requires trimming'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def trim(args):
    """
    %prog trim fastqfile

    CA6.1 cannot handle long fastq seq and requires trimming.
    """
    p = OptionParser(trim.__doc__)
    p.add_option("-l", dest="length", default=103, type="int",
            help="keep first N bases [default: %default]")

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    fastqfile, = args
    length = opts.length
    trimmedfastqfile = fastqfile.rsplit(".", 1)[0] + ".{0}.fastq".format(length)

    cmd = "fastx_trimmer -f 1 -l {0} -Q33 ".format(length)
    cmd += "< {0} > {1}".format(fastqfile, trimmedfastqfile)
    sh(cmd)
 

def make_qual(fastafile, defaultqual=20):
    """
    Make a qualfile with default qual value if not available
    """
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
    assumed as adjacent sequence records (i.e. /1, /2, /1, /2 ...).
    """
    p = OptionParser(fasta.__doc__)
    add_size_option(p)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    fastafile, = args
    libname = op.basename(fastafile).split(".")[0]
    frgfile = libname + ".frg"

    mated = (opts.size != 0)
    mean, sv = get_mean_sv(opts.size)

    qualfile = make_qual(fastafile)
    if mated:
        matefile = make_matepairs(fastafile)

    cmd = "convert-fasta-to-v2.pl -l {0} -s {1} -q {2}".format(libname,
            fastafile, qualfile)
    if mated:
        cmd += " -mean {0} -stddev {1} -m {2}".format(mean, sv, matefile)

    cmd += " > {0}".format(frgfile)
    sh(cmd)


def sff(args):
    """
    %prog sff sffiles

    Convert reads formatted as 454 SFF file, and convert to CA frg file. 
    """
    p = OptionParser(sff.__doc__)
    p.add_option("-o", dest="output", default="out",
            help="output frg filename")
    add_size_option(p)

    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(p.print_help())

    sffiles = args
    libname = opts.output

    mated = (opts.size != 0)
    mean, sv = get_mean_sv(opts.size)

    cmd = "sffToCA -libraryname {0} -output {0}".format(libname)
    cmd += " -clear 454 -trim chop"
    if mated:
        cmd += " -linker titanium -insertsize {0} {1}".format(mean, sv)

    cmd += " " + " ".join(sffiles)

    sh(cmd)
    

def fastq(args):
    """
    %prog fastq fastqfile

    Convert reads formatted as FASTQ file, and convert to CA frg file. 
    """
    p = OptionParser(fastq.__doc__)
    p.add_option("--sanger", dest="sanger", default=False, action="store_true",
            help="Are the qv sanger encodings [default: %default]?")
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

    cmd = "~/bin/Linux-amd64/bin/fastqToCA -libraryname {0} -fastq {1}".format(libname, fastqfile)
    if opts.sanger:
        cmd += "-type sanger "
    if mated:
        assert len(args)==2, "you need two fastq file for mated library"
        cmd += ",{0} -insertsize {1} {2}".format(fastqfile2, mean, sv)

    cmd += " > {0}".format(frgfile)
    sh(cmd)


if __name__ == '__main__':
    main()
