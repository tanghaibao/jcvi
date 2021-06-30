#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run bowtie2 command and skips the manual run of naming intermediate output
files. Bowtie2 help:

<http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>
"""
import sys
import logging

from jcvi.formats.base import BaseFile
from jcvi.utils.cbook import percentage
from jcvi.formats.sam import output_bam, get_prefix, get_samfile
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh, get_abs_path


first_tag = lambda fp: next(fp).split()[0]


class BowtieLogFile(BaseFile):
    """
    Simple file that contains mapping rate:

    100000 reads; of these:
      100000 (100.00%) were unpaired; of these:
        88453 (88.45%) aligned 0 times
        9772 (9.77%) aligned exactly 1 time
        1775 (1.77%) aligned >1 times
    11.55% overall alignment rate
    """

    def __init__(self, filename):

        super(BowtieLogFile, self).__init__(filename)
        fp = open(filename)
        self.total = int(first_tag(fp))
        self.unpaired = int(first_tag(fp))
        self.unmapped = int(first_tag(fp))
        self.unique = int(first_tag(fp))
        self.multiple = int(first_tag(fp))
        self.mapped = self.unique + self.multiple
        self.rate = float(first_tag(fp).rstrip("%"))
        fp.close()

    def __str__(self):
        return "Total mapped: {0}".format(percentage(self.mapped, self.total))

    __repr__ = __str__


def main():

    actions = (
        ("index", "wraps bowtie2-build"),
        ("align", "wraps bowtie2"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def check_index(dbfile):
    dbfile = get_abs_path(dbfile)
    safile = dbfile + ".1.bt2"
    if need_update(dbfile, safile):
        cmd = "bowtie2-build {0} {0}".format(dbfile)
        sh(cmd)
    else:
        logging.error("`{0}` exists. `bowtie2-build` already run.".format(safile))

    return dbfile


def index(args):
    """
    %prog index database.fasta

    Wrapper for `bowtie2-build`. Same interface.
    """
    p = OptionParser(index.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (dbfile,) = args
    check_index(dbfile)


def align(args):
    """
    %prog align database.fasta read1.fq [read2.fq]

    Wrapper for `bowtie2` single-end or paired-end, depending on the number of args.
    """
    from jcvi.formats.fastq import guessoffset

    p = OptionParser(align.__doc__)
    p.set_firstN(firstN=0)
    p.add_option(
        "--full",
        default=False,
        action="store_true",
        help="Enforce end-to-end alignment [default: local]",
    )
    p.add_option(
        "--reorder",
        default=False,
        action="store_true",
        help="Keep the input read order",
    )
    p.add_option(
        "--null",
        default=False,
        action="store_true",
        help="Do not write to SAM/BAM output",
    )
    p.add_option(
        "--fasta", default=False, action="store_true", help="Query reads are FASTA"
    )
    p.set_cutoff(cutoff=800)
    p.set_mateorientation(mateorientation="+-")
    p.set_sam_options(bowtie=True)

    opts, args = p.parse_args(args)
    extra = opts.extra
    mo = opts.mateorientation
    if mo == "+-":
        extra += ""
    elif mo == "-+":
        extra += "--rf"
    else:
        extra += "--ff"

    PE = True
    if len(args) == 2:
        logging.debug("Single-end alignment")
        PE = False
    elif len(args) == 3:
        logging.debug("Paired-end alignment")
    else:
        sys.exit(not p.print_help())

    firstN = opts.firstN
    mapped = opts.mapped
    unmapped = opts.unmapped
    fasta = opts.fasta
    gl = "--end-to-end" if opts.full else "--local"

    dbfile, readfile = args[0:2]
    dbfile = check_index(dbfile)
    prefix = get_prefix(readfile, dbfile)
    samfile, mapped, unmapped = get_samfile(
        readfile, dbfile, bowtie=True, mapped=mapped, unmapped=unmapped, bam=opts.bam
    )
    logfile = prefix + ".log"
    if not fasta:
        offset = guessoffset([readfile])

    if not need_update(dbfile, samfile):
        logging.error("`{0}` exists. `bowtie2` already run.".format(samfile))
        return samfile, logfile

    cmd = "bowtie2 -x {0}".format(dbfile)
    if PE:
        r1, r2 = args[1:3]
        cmd += " -1 {0} -2 {1}".format(r1, r2)
        cmd += " --maxins {0}".format(opts.cutoff)
        mtag, utag = "--al-conc", "--un-conc"
    else:
        cmd += " -U {0}".format(readfile)
        mtag, utag = "--al", "--un"

    if mapped:
        cmd += " {0} {1}".format(mtag, mapped)
    if unmapped:
        cmd += " {0} {1}".format(utag, unmapped)

    if firstN:
        cmd += " --upto {0}".format(firstN)
    cmd += " -p {0}".format(opts.cpus)
    if fasta:
        cmd += " -f"
    else:
        cmd += " --phred{0}".format(offset)
    cmd += " {0}".format(gl)
    if opts.reorder:
        cmd += " --reorder"

    cmd += " {0}".format(extra)
    # Finally the log
    cmd += " 2> {0}".format(logfile)

    if opts.null:
        samfile = "/dev/null"

    cmd = output_bam(cmd, samfile)
    sh(cmd)
    print(open(logfile).read(), file=sys.stderr)

    return samfile, logfile


if __name__ == "__main__":
    main()
