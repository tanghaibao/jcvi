#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Trinity assembly of RNAseq reads. Contains de novo (DN) method and genome-guided
(GG) method.

DN: http://trinityrnaseq.sourceforge.net/
GG: http://trinityrnaseq.sourceforge.net/genome_guided_trinity.html
"""

import os.path as op
import os
import sys

from jcvi.formats.base import FileMerger, write_file
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, glob


def main():

    actions = (
        ('dn', 'run trinity-dn on a folder of reads'),
        ('gg', 'run trinity-gg on a folder of reads or BAM file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def dn(args):
    """
    %prog dn folder

    Run Trinity-DN on a folder of reads. When paired-end (--paired) mode is on,
    filenames will be scanned based on whether they contain "_1_" and "_2_".
    """
    p = OptionParser(dn.__doc__)
    p.add_option("--paired", default=False, action="store_true",
                 help="Paired-end mode [default: %default]")
    p.set_trinity_opts()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    folder, = args
    paired = opts.paired
    thome = opts.trinity_home
    tfolder = folder + "_DN"

    cwd = os.getcwd()
    mkdir(tfolder)
    os.chdir(tfolder)

    flist = glob("../" + folder + "/*")
    if paired:
        f1 = [x for x in flist if "_1_" in x or ".1." in x]
        f2 = [x for x in flist if "_2_" in x or ".2." in x]
        assert len(f1) == len(f2)
        r1, r2 = "left.fastq", "right.fastq"
        reads = ((f1, r1), (f2, r2))
    else:
        r = "single.fastq"
        reads = ((flist, r), )

    for fl, r in reads:
        fm = FileMerger(fl, r)
        fm.merge(checkexists=True)

    cmd = op.join(thome, "Trinity")
    cmd += " --seqType fq --JM {0} --CPU {1}".format(opts.JM, opts.cpus)
    if paired:
        cmd += " --left {0} --right {1}".format(reads[0][-1], reads[1][-1])
    else:
        cmd += " --single {0}".format(reads[0][-1])

    runfile = "run.sh"
    write_file(runfile, cmd)
    os.chdir(cwd)


def gg(args):
    """
    %prog gg [--options] folder|bamfile genome.fasta

    Run Trinity-GG on a folder of reads or a coord-sorted BAM file.  When paired-end
    (--paired) mode is on, filenames will be scanned based on whether they contain
    the patterns ("_1_" and "_2_") or (".1." and ".2.") or ("_1." and "_2.")
    """
    p = OptionParser(gg.__doc__)
    p.add_option("--paired", default=False, action="store_true",
                 help="Paired-end mode [default: %default]")
    p.set_trinity_opts(gg=True)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    inparam, genome, = args
    paired = opts.paired
    thome = opts.trinity_home

    pf = inparam
    if not op.isdir(inparam):
        pf = inparam.split(".")[0]
    tfolder = "{0}_DN".format(pf)

    cwd = os.getcwd()
    mkdir(tfolder)
    os.chdir(tfolder)

    if not op.isdir(inparam):
        flist = glob("../" + inparam + "/*")
        if paired:
            f1 = [x for x in flist if "_1_" in x or ".1." in x or "_1." in x]
            f2 = [x for x in flist if "_2_" in x or ".2." in x or "_2." in x]
            assert len(f1) == len(f2)

    cmd = op.join(thome, "Trinity")
    cmd += " --seqType fq --JM {0} --CPU {1}".format(opts.JM, opts.cpus)
    cmd += " --genome {0} --genome_guided_max_intron {1}".format(genome, opts.max_intron)
    if flist:
        if paired:
            cmd += " --left {0} --right {1}".format(" ".join(f1), " ".join(f2))
        else:
            cmd += " --single {0}".format(" ".join(flist))
    else:
        cmd += " --genome_guided_use_bam {0}".format(inparam)

    runfile = "run.sh"
    write_file(runfile, cmd)
    os.chdir(cwd)


if __name__ == '__main__':
    main()
