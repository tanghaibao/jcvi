#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run bwa command and skips the manual run of naming intermediate output files
The whole pipeline is following bwa documentation at
<http://bio-bwa.sf.net/bwa.shtml>
"""

import sys
import logging
import os.path as op

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, set_grid, set_params, need_update, \
                sh, debug
debug()


def main():

    actions = (
        ('index', 'wraps bwa index'),
        ('aln', 'wraps bwa aln'),
        ('samse', 'wraps bwa samse'),
        ('sampe', 'wraps bwa sampe'),
        ('bwasw', 'wraps bwa bwasw'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def check_index(dbfile, grid=False):
    safile = dbfile + ".sa"
    if need_update(dbfile, safile):
        cmd = "bwa index -a bwtsw {0}".format(dbfile)
        sh(cmd, grid=grid)
    else:
        logging.error("`{0}` exists. `bwa index` already run.".format(safile))

    return safile


def check_aln(dbfile, readfile, grid=False, cpus=32):
    from jcvi.formats.fastq import guessoffset

    saifile = readfile.rsplit(".", 1)[0] + ".sai"
    if need_update((dbfile, readfile), saifile):
        offset = guessoffset([readfile])
        cmd = "bwa aln -t {0}".format(cpus)
        if offset == 64:
            cmd += " -I"

        cmd += " {0} {1}".format(dbfile, readfile)
        sh(cmd, grid=grid, outfile=saifile)
    else:
        logging.error("`{0}` exists. `bwa aln` already run.".format(saifile))

    return saifile


def index(args):
    """
    %prog index database.fasta

    Wrapper for `bwa index`. Same interface, only adds grid submission.
    """
    p = OptionParser(index.__doc__)
    set_params(p)
    set_grid(p)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    extra = opts.extra
    grid = opts.grid

    dbfile, = args
    safile = check_index(dbfile, grid=grid)


def aln(args):
    """
    %prog aln database.fasta *.fastq

    Wrapper for `bwa aln` except this will run over a set of files.
    """
    p = OptionParser(aln.__doc__)
    p.add_option("--cpus", default=32,
                 help="Number of cpus to use [default: %default]")
    set_params(p)
    set_grid(p)

    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(p.print_help())

    extra = opts.extra
    grid = opts.grid

    dbfile, readfiles = args[0], args[1:]
    safile = check_index(dbfile, grid=grid)
    for readfile in readfiles:
        saifile = check_aln(dbfile, readfile, grid=grid, cpus=opts.cpus)


def samse(args):
    """
    %prog samse database.fasta short_read.fastq

    Wrapper for `bwa samse`. Output will be short_read.sam.
    """
    p = OptionParser(samse.__doc__)
    p.add_option("--bam", default=False, action="store_true",
                 help="write to bam file [default: %default]")
    p.add_option("--cpus", default=32,
                 help="Number of cpus to use [default: %default]")
    set_params(p)
    set_grid(p)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    extra = opts.extra
    grid = opts.grid

    dbfile, readfile = args
    safile = check_index(dbfile, grid=grid)
    saifile = check_aln(dbfile, readfile, grid=grid, cpus=opts.cpus)

    prefix = readfile.rsplit(".", 1)[0]
    samfile = (prefix + ".bam") if opts.bam else (prefix + ".sam")
    if not need_update((safile, saifile), samfile):
        logging.error("`{0}` exists. `bwa samse` already run.".format(samfile))
        return

    cmd = "bwa samse {0} {1} {2} ".format(dbfile, saifile, readfile)
    cmd += "{0}".format(extra)
    if opts.bam:
        cmd += " | samtools view -bS -F 4 - "
    sh(cmd, grid=grid, outfile=samfile, threaded=opts.cpus)


def sampe(args):
    """
    %prog sampe database.fasta read1.fq read2.fq

    Wrapper for `bwa sampe`. Output will be read1.sam.
    """
    p = OptionParser(sampe.__doc__)
    p.add_option("--bam", default=False, action="store_true",
                 help="write to bam file [default: %default]")
    p.add_option("--cpus", default=32,
                 help="Number of cpus to use [default: %default]")
    set_params(p)
    set_grid(p)

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(p.print_help())

    extra = opts.extra
    grid = opts.grid

    dbfile, read1file, read2file = args
    safile = check_index(dbfile, grid=grid)
    sai1file = check_aln(dbfile, read1file, grid=grid, cpus=opts.cpus)
    sai2file = check_aln(dbfile, read2file, grid=grid, cpus=opts.cpus)

    prefix = read1file.rsplit(".", 1)[0]
    samfile = (prefix + ".bam") if opts.bam else (prefix + ".sam")
    if not need_update((safile, sai1file, sai2file), samfile):
        logging.error("`{0}` exists. `bwa samse` already run.".format(samfile))
        return

    cmd = "bwa sampe {0} {1} {2} {3} {4} ".format(dbfile, sai1file, sai2file,
            read1file, read2file)
    cmd += "{0}".format(extra)
    if opts.bam:
        cmd += " | samtools view -bS -F 4 - "
    sh(cmd, grid=grid, outfile=samfile)


def bwasw(args):
    """
    %prog bwasw database.fasta long_read.fastq

    Wrapper for `bwa bwasw`. Output will be long_read.sam.
    """
    p = OptionParser(bwasw.__doc__)
    set_params(p)
    set_grid(p)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    extra = opts.extra
    grid = opts.grid

    dbfile, readfile = args
    safile = check_index(dbfile, grid=grid)
    saifile = check_aln(dbfile, readfile, grid=grid)

    samfile = readfile.rsplit(".", 1)[0] + ".sam"
    if not need_update((safile, saifile), samfile):
        logging.error("`{0}` exists. `bwa bwasw` already run.".format(samfile))
        return

    cmd = "bwa bwasw -t 32 {0} {1} ".format(dbfile, readfile)
    cmd += "{0}".format(extra)
    sh(cmd, grid=grid, outfile=samfile)


if __name__ == '__main__':
    main()
