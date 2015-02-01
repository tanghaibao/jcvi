#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run bwa command and skips the manual run of naming intermediate output files
The whole pipeline is following bwa documentation at
<http://bio-bwa.sf.net/bwa.shtml>
"""

import sys
import logging

from jcvi.formats.sam import output_bam, get_samfile, mapped
from jcvi.formats.base import FileShredder
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh, \
            get_abs_path


def main():

    actions = (
        ('index', 'wraps bwa index'),
        ('align', 'wraps bwa aln|mem|bwasw'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def check_index(dbfile):
    dbfile = get_abs_path(dbfile)
    safile = dbfile + ".sa"
    if need_update(dbfile, safile):
        cmd = "bwa index {0}".format(dbfile)
        sh(cmd)
    else:
        logging.error("`{0}` exists. `bwa index` already run.".format(safile))

    return dbfile


def check_aln(dbfile, readfile, cpus=32):
    from jcvi.formats.fastq import guessoffset

    saifile = readfile.rsplit(".", 1)[0] + ".sai"
    if need_update((dbfile, readfile), saifile):
        offset = guessoffset([readfile])
        cmd = "bwa aln " + " ".join((dbfile, readfile))
        cmd += " -t {0}".format(cpus)
        if offset == 64:
            cmd += " -I"
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
    check_index(dbfile)


def align(args):
    """
    %prog align database.fasta read1.fq [read2.fq]

    Wrapper for three modes of BWA - mem (default), aln, bwasw (long reads).
    """
    valid_modes = ("bwasw", "aln", "mem")
    p = OptionParser(align.__doc__)
    p.add_option("--mode", default="mem", choices=valid_modes, help="BWA mode")
    p.add_option("--readtype", choices=("pacbio", "pbread", "ont2d", "intractg"),
                 help="Read type in bwa-mem")
    p.set_cutoff(cutoff=800)
    p.set_sam_options()

    opts, args = p.parse_args(args)
    mode = opts.mode
    nargs = len(args)

    if nargs not in (2, 3):
        sys.exit(not p.print_help())

    tag = "bwa-{0}: ".format(mode)
    c = mem
    if nargs == 2:
        tag += "Single-end alignment"
        if mode == "bwasw":
            c = bwasw
        elif mode == "aln":
            c = samse
    else:
        assert mode != "bwasw", "Cannot use --bwasw with paired-end mode"
        tag += "Paired-end alignment"
        if mode == "aln":
            c = sampe

    logging.debug(tag)
    cmd, samfile = c(args, opts)
    if cmd:
        cmd = output_bam(cmd, samfile)

    bam = opts.bam
    unmapped = opts.unmapped

    sh(cmd)
    if unmapped:
        dbfile, readfile = args[:2]
        mopts = [samfile, "--unmapped"]
        if not bam:
            mopts += ["--sam"]
        mapped(mopts)
        FileShredder([samfile])

    return samfile, None


def samse(args, opts):
    """
    %prog samse database.fasta short_read.fastq

    Wrapper for `bwa samse`. Output will be short_read.sam.
    """
    dbfile, readfile = args
    dbfile = check_index(dbfile)
    saifile = check_aln(dbfile, readfile, cpus=opts.cpus)

    samfile, _, unmapped = get_samfile(readfile, dbfile,
                                       bam=opts.bam, unmapped=opts.unmapped)
    if not need_update((dbfile, saifile), samfile):
        logging.error("`{0}` exists. `bwa samse` already run.".format(samfile))
        return "", samfile

    cmd = "bwa samse {0} {1} {2}".format(dbfile, saifile, readfile)
    cmd += " " + opts.extra
    if opts.uniq:
        cmd += " -n 1"

    return cmd, samfile


def sampe(args, opts):
    """
    %prog sampe database.fasta read1.fq read2.fq

    Wrapper for `bwa sampe`. Output will be read1.sam.
    """
    dbfile, read1file, read2file = args
    dbfile = check_index(dbfile)
    sai1file = check_aln(dbfile, read1file, cpus=opts.cpus)
    sai2file = check_aln(dbfile, read2file, cpus=opts.cpus)

    samfile, _, unmapped = get_samfile(read1file, dbfile,
                                       bam=opts.bam, unmapped=opts.unmapped)
    if not need_update((dbfile, sai1file, sai2file), samfile):
        logging.error("`{0}` exists. `bwa samse` already run.".format(samfile))
        return "", samfile

    cmd = "bwa sampe " + " ".join((dbfile, sai1file, sai2file, \
                                   read1file, read2file))
    cmd += " " + opts.extra
    if opts.cutoff:
        cmd += " -a {0}".format(opts.cutoff)
    if opts.uniq:
        cmd += " -n 1"

    return cmd, samfile


def mem(args, opts):
    """
    %prog mem database.fasta read1.fq [read2.fq]

    Wrapper for `bwa mem`. Output will be read1.sam.
    """
    dbfile, read1file = args[:2]

    dbfile = check_index(dbfile)
    args[0] = dbfile
    samfile, _, unmapped = get_samfile(read1file, dbfile,
                                       bam=opts.bam, unmapped=opts.unmapped)
    if not need_update(read1file, samfile):
        logging.error("`{0}` exists. `bwa mem` already run.".format(samfile))
        return "", samfile

    cmd = "bwa mem " + " ".join(args)
    cmd += " -M -t {0}".format(opts.cpus)
    if opts.readtype:
        cmd += " -x {0}".format(opts.readtype)
    cmd += " " + opts.extra
    return cmd, samfile


def bwasw(args, opts):
    """
    %prog bwasw database.fasta long_read.fastq

    Wrapper for `bwa bwasw`. Output will be long_read.sam.
    """
    dbfile, readfile = args
    dbfile = check_index(dbfile)

    samfile, _, unmapped = get_samfile(readfile, dbfile,
                                       bam=opts.bam, unmapped=opts.unmapped)
    if not need_update(dbfile, samfile):
        logging.error("`{0}` exists. `bwa bwasw` already run.".format(samfile))
        return "", samfile

    cmd = "bwa bwasw " + " ".join(args)
    cmd += " -t {0}".format(opts.cpus)
    cmd += " " + opts.extra
    return cmd, samfile


if __name__ == '__main__':
    main()
