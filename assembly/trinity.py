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
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, iglob


def main():

    actions = (
        ('prepare', 'prepare shell script to run trinity-dn/gg on a folder of reads'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def prepare(args):
    """
    %prog prepare [--options] folder [genome.fasta]

    Run Trinity on a folder of reads.  When paired-end (--paired) mode is on,
    filenames will be scanned based on whether they contain the patterns
    ("_1_" and "_2_") or (".1." and ".2.") or ("_1." and "_2.").

    By default, prepare script for DN

    If genome.fasta is provided, prepare script for GG-Trinity.
    If coord-sorted BAM is provided, then it will use it as starting point.

    Since GG-Trinity jobs are partitioned DN-Trinity jobs run on relatively small
    regions, lesser amount of CPU can be specified for each DN job using `--gg_cpu`
    In such cases, the `--cpu` should be set to a larger value to help speedup
    upstream steps such as GSNAP read mapping or coordinate sorting of BAM files.

    Newer versions of trinity can take multiple fastq files as input.
    If "--merge" is specified, the fastq files are merged together before assembling
    """
    p = OptionParser(prepare.__doc__)
    p.add_option("--paired", default=False, action="store_true",
                 help="Paired-end mode [default: %default]")
    p.add_option("--merge", default=False, action="store_true",
                 help="Merge individual input fastq's into left/right/single" + \
                      " file(s) [default: %default]")
    p.set_trinity_opts()
    p.set_grid()
    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    inparam, = args[:1]
    genome = args[1] if len(args) == 2 else None
    method = "GG" if genome is not None else "DN"

    paired = opts.paired
    merge = opts.merge
    thome = opts.trinity_home
    use_bam = opts.use_bam
    gg_cpu = opts.gg_cpu

    pf = inparam.split(".")[0]
    tfolder = "{0}_{1}".format(pf, method)

    cwd = os.getcwd()
    mkdir(tfolder)
    os.chdir(tfolder)

    flist = iglob("../" + inparam, "*.fq", "*.fastq", "*.fq.gz", "*.fastq.gz")
    if paired:
        f1 = [x for x in flist if "_1_" in x or ".1." in x or "_1." in x]
        f2 = [x for x in flist if "_2_" in x or ".2." in x or "_2." in x]
        assert len(f1) == len(f2)
        if merge:
            r1, r2 = "left.fastq", "right.fastq"
            reads = ((f1, r1), (f2, r2))
    else:
        if merge:
            r = "single.fastq"
            reads = ((flist, r), )

    if merge:
        for fl, r in reads:
            fm = FileMerger(fl, r)
            fm.merge(checkexists=True)

    cmd = op.join(thome, "Trinity")
    cmd += " --seqType fq --JM {0} --CPU {1}".format(opts.JM, opts.cpus)
    cmd += " --min_contig_length {0}".format(opts.min_contig_length)
    if opts.bflyGCThreads:
        cmd += " --bflyGCThreads {0}".format(opts.bflyGCThreads)

    if method == "GG":
        cmd += " --genome {0} --genome_guided_max_intron {1}".format(genome, opts.max_intron)
        if use_bam:
            cmd += " --genome_guided_use_bam {0}".format(use_bam)
        if gg_cpu:
            cmd += " --genome_guided_CPU {0}".format(gg_cpu)
    if opts.grid and opts.grid_conf_file:
        cmd += " --grid_conf_file={0}".format(opts.grid_conf_file)

    if paired:
        if merge:
            cmd += " --left {0} --right {1}".format(reads[0][-1], reads[1][-1])
        else:
            for lf, rf in zip(f1, f2):
                cmd += " --left {0}".format(lf)
                cmd += " --right {0}".format(rf)
    else:
        if merge:
             cmd += " --single {0}".format(reads[0][-1])
        else:
            for f in flist:
                cmd += " --single {0}".format(f)
    if opts.extra:
        cmd += " {0}".format(opts.extra)

    runfile = "run.sh"
    write_file(runfile, cmd)
    os.chdir(cwd)


if __name__ == '__main__':
    main()
