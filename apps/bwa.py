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

from jcvi.formats.sam import output_bam, get_samfile, mapped
from jcvi.formats.base import FileShredder
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, \
                sh, debug
debug()


def main():

    actions = (
        ('index', 'wraps bwa index'),
        ('align', 'wraps bwa samse or sampe'),
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
    p.set_params()
    p.set_grid()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    extra = opts.extra
    grid = opts.grid

    dbfile, = args
    safile = check_index(dbfile, grid=grid)


def align(args):
    """
    %prog align database.fasta read1.fq [read2.fq]

    Wrapper for `bwa samse` or `bwa sampe`, depending on the number of args.
    """
    p = OptionParser(align.__doc__)
    p.add_option("--bwasw", default=False, action="store_true",
                 help="Run bwasw mode [default: %default]")
    p.set_sam_options()

    opts, args = p.parse_args(args)
    bsw = opts.bwasw

    if len(args) not in (2, 3):
        sys.exit(not p.print_help())

    if len(args) == 2:
        mode = "Single-end alignment"
        if bsw:
            mode += " (long reads)"
            c = bwasw
        else:
            c = samse
    else:
        assert not bsw, "Cannot use --bwasw with paired-end mode"
        mode = "Paired-end alignment"
        c = sampe

    logging.debug(mode)
    cmd = c(args, opts)

    bam = opts.bam
    grid = opts.grid
    unmapped = opts.unmapped

    sh(cmd, grid=grid, threaded=opts.cpus)
    if grid:
        return

    if unmapped:
        dbfile, readfile = args[:2]
        samfile, _, unmapped = get_samfile(readfile, dbfile,
                                           bam=bam, unmapped=unmapped)
        mopts = [samfile, "--unmapped"]
        if bam:
            mopts += ["--bam"]
        mapped(mopts)
        FileShredder([samfile])


def samse(args, opts):
    """
    %prog samse database.fasta short_read.fastq

    Wrapper for `bwa samse`. Output will be short_read.sam.
    """
    extra = opts.extra
    grid = opts.grid

    dbfile, readfile = args
    safile = check_index(dbfile, grid=grid)
    saifile = check_aln(dbfile, readfile, grid=grid, cpus=opts.cpus)

    samfile, _, unmapped = get_samfile(readfile, dbfile,
                                       bam=opts.bam, unmapped=opts.unmapped)
    if not need_update((safile, saifile), samfile):
        logging.error("`{0}` exists. `bwa samse` already run.".format(samfile))
        return

    cmd = "bwa samse {0} {1} {2}".format(dbfile, saifile, readfile)
    cmd += " {0}".format(extra)
    if opts.uniq:
        cmd += " -n 1"

    return output_bam(cmd, samfile)


def sampe(args, opts):
    """
    %prog sampe database.fasta read1.fq read2.fq

    Wrapper for `bwa sampe`. Output will be read1.sam.
    """
    extra = opts.extra
    grid = opts.grid

    dbfile, read1file, read2file = args
    safile = check_index(dbfile, grid=grid)
    sai1file = check_aln(dbfile, read1file, grid=grid, cpus=opts.cpus)
    sai2file = check_aln(dbfile, read2file, grid=grid, cpus=opts.cpus)

    samfile, _, unmapped = get_samfile(read1file, dbfile,
                                       bam=opts.bam, unmapped=opts.unmapped)
    if not need_update((safile, sai1file, sai2file), samfile):
        logging.error("`{0}` exists. `bwa samse` already run.".format(samfile))
        return

    cmd = "bwa sampe {0} {1} {2} {3} {4}".format(dbfile, sai1file, sai2file,
            read1file, read2file)
    cmd += " {0}".format(extra)
    if opts.uniq:
        cmd += " -n 1"

    return output_bam(cmd, samfile)


def bwasw(args, opts):
    """
    %prog bwasw database.fasta long_read.fastq

    Wrapper for `bwa bwasw`. Output will be long_read.sam.
    """
    grid = opts.grid
    cpus = opts.cpus
    bam = opts.bam
    unmapped = opts.unmapped

    dbfile, readfile = args
    safile = check_index(dbfile, grid=grid)

    samfile, _, unmapped = get_samfile(readfile, dbfile,
                                       bam=bam, unmapped=unmapped)
    if not need_update(safile, samfile):
        logging.error("`{0}` exists. `bwa bwasw` already run.".format(samfile))
        return

    cmd = "bwa bwasw -t {0} {1} {2}".format(cpus, dbfile, readfile)
    cmd += "{0}".format(opts.extra)
    return output_bam(cmd, samfile)


if __name__ == '__main__':
    main()
