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
from jcvi.apps.softlink import get_abs_path
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


def check_index(dbfile):
    safile = dbfile + ".sa"
    if need_update(dbfile, safile):
        cmd = "bwa index -a bwtsw {0}".format(dbfile)
        sh(cmd)
    else:
        logging.error("`{0}` exists. `bwa index` already run.".format(safile))

    return safile


def check_aln(dbfile, readfile, cpus=32):
    from jcvi.formats.fastq import guessoffset

    saifile = readfile.rsplit(".", 1)[0] + ".sai"
    if need_update((dbfile, readfile), saifile):
        offset = guessoffset([readfile])
        cmd = "bwa aln -t {0}".format(cpus)
        if offset == 64:
            cmd += " -I"

        cmd += " {0} {1}".format(dbfile, readfile)
        sh(cmd, outfile=saifile)
    else:
        logging.error("`{0}` exists. `bwa aln` already run.".format(saifile))

    return saifile


def index(args):
    """
    %prog index database.fasta

    Wrapper for `bwa index`. Same interface.
    """
    p = OptionParser(index.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    dbfile, = args
    safile = check_index(dbfile)


def align(args):
    """
    %prog align database.fasta read1.fq [read2.fq]

    Wrapper for `bwa samse` or `bwa sampe`, depending on the number of args.
    """
    p = OptionParser(align.__doc__)
    p.add_option("--bwasw", default=False, action="store_true",
                 help="Run bwasw mode [default: %default]")
    p.set_cutoff(cutoff=800)
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

    args[0] = get_abs_path(args[0])
    logging.debug(mode)
    cmd, samfile = c(args, opts)
    if cmd:
        cmd = output_bam(cmd, samfile)

    bam = opts.bam
    unmapped = opts.unmapped

    sh(cmd, threaded=opts.cpus)
    if unmapped:
        dbfile, readfile = args[:2]
        mopts = [samfile, "--unmapped"]
        if bam:
            mopts += ["--bam"]
        mapped(mopts)
        FileShredder([samfile])

    return samfile, None


def samse(args, opts):
    """
    %prog samse database.fasta short_read.fastq

    Wrapper for `bwa samse`. Output will be short_read.sam.
    """
    dbfile, readfile = args
    safile = check_index(dbfile)
    saifile = check_aln(dbfile, readfile, cpus=opts.cpus)

    samfile, _, unmapped = get_samfile(readfile, dbfile,
                                       bam=opts.bam, unmapped=opts.unmapped)
    if not need_update((safile, saifile), samfile):
        logging.error("`{0}` exists. `bwa samse` already run.".format(samfile))
        return "", samfile

    cmd = "bwa samse {0} {1} {2}".format(dbfile, saifile, readfile)
    cmd += " {0}".format(opts.extra)
    if opts.uniq:
        cmd += " -n 1"

    return cmd, samfile


def sampe(args, opts):
    """
    %prog sampe database.fasta read1.fq read2.fq

    Wrapper for `bwa sampe`. Output will be read1.sam.
    """
    dbfile, read1file, read2file = args
    safile = check_index(dbfile)
    sai1file = check_aln(dbfile, read1file, cpus=opts.cpus)
    sai2file = check_aln(dbfile, read2file, cpus=opts.cpus)

    samfile, _, unmapped = get_samfile(read1file, dbfile,
                                       bam=opts.bam, unmapped=opts.unmapped)
    if not need_update((safile, sai1file, sai2file), samfile):
        logging.error("`{0}` exists. `bwa samse` already run.".format(samfile))
        return

    cmd = "bwa sampe {0} {1} {2} {3} {4}".format(dbfile, sai1file, sai2file,
            read1file, read2file)
    if opts.cutoff:
        cmd += " -a {0}".format(opts.cutoff)
    cmd += " {0}".format(opts.extra)
    if opts.uniq:
        cmd += " -n 1"

    return cmd, samfile


def bwasw(args, opts):
    """
    %prog bwasw database.fasta long_read.fastq

    Wrapper for `bwa bwasw`. Output will be long_read.sam.
    """
    cpus = opts.cpus
    bam = opts.bam
    unmapped = opts.unmapped

    dbfile, readfile = args
    safile = check_index(dbfile)

    samfile, _, unmapped = get_samfile(readfile, dbfile,
                                       bam=bam, unmapped=unmapped)
    if not need_update(safile, samfile):
        logging.error("`{0}` exists. `bwa bwasw` already run.".format(samfile))
        return

    cmd = "bwa bwasw -t {0} {1} {2}".format(cpus, dbfile, readfile)
    cmd += "{0}".format(opts.extra)
    return cmd, samfile


if __name__ == '__main__':
    main()
