#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper script for some programs in clc-ngs-cell
"""

import sys
import os.path as op

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug, set_grid, set_params, sh
debug()


def main():

    actions = (
        ('map', 'map the reads to the reference'),
        ('trim', 'wrapper around clc quality_trim'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def map(args):
    """
    %prog map reference fastqfiles

    Use `clc_ref_assemble` to map the read files to a reference. Use a non-zero
    -s option to turn on paired end mode.
    """
    p = OptionParser(map.__doc__)
    p.add_option("-o", dest="outfile", default=None,
            help="Output prefix.cas file [default: %default]")
    p.add_option("-s", dest="size", default=0, type="int",
            help="Use paired end mapping with insert [default: %default]")
    p.add_option("--short", default=False, action="store_true",
            help="Use `clc_ref_assemble_short` as the mapper [default: %default]")
    p.add_option("--orientations", default="fb",
            help="The reads have the orientations [default: %default]")
    set_params(p)
    set_grid(p)

    opts, args = p.parse_args(args)
    if len(args) < 2:
        sys.exit(not p.print_help())

    license = "license.properties"
    if not op.exists(license):
        sh("cp ~/{0} .".format(license))

    ref = args[0]
    assert op.exists(ref)
    fastqfiles = args[1:]
    size = opts.size
    orientations = opts.orientations
    assert orientations in ("fb", "bf", "ff", "bb")

    cmd = "clc_ref_assemble_short" if opts.short else "clc_ref_assemble_long"
    readprefix = op.basename(fastqfiles[0]).split(".", 1)[0]
    refprefix = op.basename(ref).split(".", 1)[0]
    outfile = opts.outfile or "{0}.{1}".format(readprefix, refprefix)
    if not outfile.endswith(".cas"):
        outfile += ".cas"

    cmd += " --cpus 16"
    cmd += " -d {0} -o {1} -q ".format(ref, outfile)
    fastqs = " ".join(fastqfiles)
    if size == 0:
        cmd += fastqs
    else:
        assert len(fastqfiles) == 2
        stddev = size / 4
        lb, ub = size - stddev, size + stddev
        cmd += " -p {0} ss {1} {2} -i {3} ".format(orientations, lb, ub, fastqs)

    if opts.extra:
        cmd += " " + opts.extra

    if not opts.short:
        cmd += " -l 0.9 -s 0.9"

    sh(cmd, grid=opts.grid)


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
    p.add_option("--offset", dest="offset", type="int", default=64,
            help="Set the ascii offset value in fastq [default: %default]")
    p.add_option("--fasta", dest="fasta", default=False, action="store_true",
            help="Output fasta sequence? [default: fastq]")
    set_grid(p)

    opts, args = p.parse_args(args)

    largs = len(args)
    if largs not in (1, 2):
        sys.exit(p.print_help())

    paired = (largs == 2)
    fastqfile1 = args[0]
    assert op.exists(fastqfile1)

    suffix = "fasta" if opts.fasta else "fastq"

    if paired:
        fastqfile2 = args[1]
        assert op.exists(fastqfile2)

    prefix = fastqfile1.split('.')[0]
    cmd = "quality_trim -c {0.cutoff} -m {0.minlength} -f {0.offset} ".format(opts)
    if paired:
        cmd += "-r -i {0} {1} ".format(fastqfile1, fastqfile2)
        cmd += "-p {0}.pairs.{1} ".format(prefix, suffix)
    else:
        cmd += "-r {0} ".format(fastqfile1)

    cmd += "-o {0}.fragments.{1}".format(prefix, suffix)
    sh(cmd, grid=opts.grid)


if __name__ == '__main__':
    main()
