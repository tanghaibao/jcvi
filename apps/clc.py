#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper script for some programs in clc-ngs-cell
"""

import sys
import os.path as op

from optparse import OptionParser

from jcvi.apps.grid import GridProcess
from jcvi.apps.base import ActionDispatcher, debug, set_grid, sh
debug()


def main():

    actions = (
        ('trim', 'wrapper around clc quality_trim'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def trim(args):
    """
    %prog trim fastqfiles

    Use `quality_trim` to trim fastq files. If there are two fastqfiles
    inputted, it is assumed as pairs of fastqs.
    """
    p = OptionParser(trim.__doc__)

    # There are many more options from `quality_trim`, but most useful twos are
    # quality cutoff (-c) and length cutoff (-m)
    p.add_option("-c", "--cutoff", dest="cutoff", type="int", default=20,
            help="Set the minimum quality for a good nucleotide. " +\
                 "[default: %default]")
    p.add_option("-m", "--minlength", dest="minlength", type="int", default=30,
            help="Set the minimum length of output reads. " +\
                 "[default: %default]")
    p.add_option("--fasta", dest="fasta", default=False, action="store_true",
            help="Output fasta sequence? [default: fastq]")
    set_grid(p, args)

    opts, args = p.parse_args(args)
    
    largs = len(args)
    if largs not in (1, 2):
        sys.exit(p.print_help())

    paired = (largs==2)
    fastqfile1 = args[0]
    assert op.exists(fastqfile1)

    suffix = "fasta" if opts.fasta else "fastq"

    if paired: 
        fastqfile2 = args[1]
        assert op.exists(fastqfile2)

    prefix = fastqfile1.split('.')[0]
    cmd = "quality_trim -c {0.cutoff} -m {0.minlength} ".format(opts)
    if paired:
        cmd += "-r -i {0} {1} ".format(fastqfile1, fastqfile2)
        cmd += "-p {0}.pairs.{1} ".format(prefix, suffix)
    else:
        cmd += "-r {0} ".format(fastqfile1)

    cmd += "-o {0}.fragments.{1}".format(prefix, suffix)

    if opts.grid:
        pr = GridProcess(cmd)
        pr.start(path=None)
    else:
        sh(cmd)


if __name__ == '__main__':
    main()
