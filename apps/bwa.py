#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run bwa command and skips the manual run of naming intermediate output files
The whole pipeline is following bwa documentation at
<http://bio-bwa.sf.net/bwa.shtml>
"""

import os.path as op
import sys
import logging

from jcvi.formats.sam import output_bam, get_samfile, mapped
from jcvi.formats.base import FileShredder
from jcvi.assembly.automaton import iter_project
from jcvi.apps.grid import MakeManager
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh, \
            get_abs_path, mkdir


def main():

    actions = (
        ('index', 'wraps bwa index'),
        ('align', 'wraps bwa aln|mem|bwasw'),
        ('batch', 'run bwa in batch mode'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def batch(args):
    """
    %proj batch database.fasta project_dir output_dir

    Run bwa in batch mode.
    """
    p = OptionParser(batch.__doc__)
    set_align_options(p)
    p.set_sam_options()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    ref_fasta, proj_dir, outdir = args
    outdir = outdir.rstrip("/")
    s3dir = None
    if outdir.startswith("s3://"):
        s3dir = outdir
        outdir = op.basename(outdir)
        mkdir(outdir)

    mm = MakeManager()
    for p, pf in iter_project(proj_dir):
        targs = [ref_fasta] + p
        cmd1, bamfile = mem(targs, opts)
        if cmd1:
            cmd1 = output_bam(cmd1, bamfile)
        nbamfile = op.join(outdir, bamfile)
        cmd2 = "mv {} {}".format(bamfile, nbamfile)
        cmds = [cmd1, cmd2]

        if s3dir:
            cmd = "aws s3 cp {} {} --sse".format(nbamfile,
                                              op.join(s3dir, bamfile))
            cmds.append(cmd)

        mm.add(p, nbamfile, cmds)

    mm.write()


def check_index(dbfile):
    dbfile = get_abs_path(dbfile)
    safile = dbfile + ".sa"
    if not op.exists(safile):
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


def set_align_options(p):
    """ Used in align() and batch()
    """
    p.add_option("--bwa", default="bwa", help="Run bwa at this path")
    p.add_option("--rg", help="Read group")
    p.add_option("--readtype",
                 choices=("pacbio", "pbread", "ont2d", "intractg"),
                 help="Read type in bwa-mem")
    p.set_cutoff(cutoff=800)


def align(args):
    """
    %prog align database.fasta read1.fq [read2.fq]

    Wrapper for three modes of BWA - mem (default), aln, bwasw (long reads).
    """
    valid_modes = ("bwasw", "aln", "mem")
    p = OptionParser(align.__doc__)
    p.add_option("--mode", default="mem", choices=valid_modes, help="BWA mode")
    set_align_options(p)
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

    cmd = "bwa sampe " + " ".join((dbfile, sai1file, sai2file,
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
    readtype = opts.readtype
    pl = readtype or "illumina"

    pf = op.basename(read1file).split(".")[0]
    rg = opts.rg or r"@RG\tID:{0}\tSM:sm\tLB:lb\tPL:{1}".format(pf, pl)
    dbfile = check_index(dbfile)
    args[0] = dbfile
    samfile, _, unmapped = get_samfile(read1file, dbfile,
                                       bam=opts.bam, unmapped=opts.unmapped)
    if not need_update(read1file, samfile):
        logging.error("`{0}` exists. `bwa mem` already run.".format(samfile))
        return "", samfile

    cmd = "{} mem".format(opts.bwa)
    '''
    -M Mark shorter split hits as secondary (for Picard compatibility).
    '''
    cmd += " -M -t {0}".format(opts.cpus)
    cmd += ' -R "{0}"'.format(rg)
    if readtype:
        cmd += " -x {0}".format(readtype)
    cmd += " " + opts.extra
    cmd += " ".join(args)

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
