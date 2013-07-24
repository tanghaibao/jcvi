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

from jcvi.formats.sam import output_bam, add_sam_options
from jcvi.apps.base import ActionDispatcher, set_grid, set_params, need_update, \
                sh, debug
debug()


def main():

    actions = (
        ('index', 'wraps bowtie2-build'),
        ('align', 'wraps bowtie2'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def check_index(dbfile, grid=False):
    safile = dbfile + ".1.bt2"
    if need_update(dbfile, safile):
        cmd = "bowtie2-build {0} {0}".format(dbfile)
        sh(cmd, grid=grid)
    else:
        logging.error("`{0}` exists. `bowtie2-build` already run.".format(safile))

    return safile


def index(args):
    """
    %prog index database.fasta

    Wrapper for `bowtie2-build`. Same interface, only adds grid submission.
    """
    p = OptionParser(index.__doc__)
    set_params(p)
    set_grid(p)

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

    Wrapper for `bowtie2` single-end or paired-end, depending on the number of args.
    """
    from jcvi.formats.fastq import guessoffset

    p = OptionParser(align.__doc__)
    add_sam_options(p)

    opts, args = p.parse_args(args)
    extra = opts.extra
    grid = opts.grid

    PE = True
    if len(args) == 2:
        logging.debug("Single-end alignment")
        PE = False
    elif len(args) == 3:
        logging.debug("Paired-end alignment")
    else:
        sys.exit(not p.print_help())

    extra = opts.extra
    grid = opts.grid

    dbfile, readfile = args[0:2]
    safile = check_index(dbfile, grid=grid)

    prefix = readfile.rsplit(".", 1)[0]
    samfile = (prefix + ".bam") if opts.bam else (prefix + ".sam")
    offset = guessoffset([readfile])

    if not need_update(safile, samfile):
        logging.error("`{0}` exists. `bowtie2` already run.".format(samfile))
        return

    cmd = "bowtie2 -x {0}".format(dbfile)
    if PE:
        r1, r2 = args[1:3]
        cmd += " -1 {0} -2 {1}".format(r1, r2)
    else:
        cmd += " -U {0}".format(readfile)
    cmd += " -p {0}".format(opts.cpus)
    cmd += " --phred{0}".format(offset)
    cmd += " {0}".format(extra)

    cmd = output_bam(cmd, bam=opts.bam)
    sh(cmd, grid=grid, outfile=samfile, threaded=opts.cpus)


if __name__ == '__main__':
    main()
