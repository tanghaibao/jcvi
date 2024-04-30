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
        ("prepare", "prepare shell script to run trinity-dn/gg on a folder of reads"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def prepare(args):
    """
    %prog prepare [--options] folder [--bam rnaseq.coordSorted.bam]

    Run Trinity on a folder of reads.  When paired-end (--paired) mode is on,
    filenames will be scanned based on whether they contain the patterns
    ("_1_" and "_2_") or (".1." and ".2.") or ("_1." and "_2.").

    By default, prepare script for DN-Trinity.

    If coord-sorted BAM is provided, prepare script for GG-Trinity, using BAM
    as starting point.

    Newer versions of trinity can take multiple fastq files as input.
    If "--merge" is specified, the fastq files are merged together before assembling
    """
    p = OptionParser(prepare.__doc__)
    p.add_option(
        "--paired",
        default=False,
        action="store_true",
        help="Paired-end mode",
    )
    p.add_option(
        "--merge",
        default=False,
        action="store_true",
        help="Merge individual input fastq's into left/right/single file(s)",
    )
    p.set_trinity_opts()
    p.set_fastq_names()
    p.set_grid()
    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    (inparam,) = args[:1]

    paired = opts.paired
    merge = opts.merge
    trinity_home = opts.trinity_home
    hpc_grid_runner_home = opts.hpcgridrunner_home

    method = "DN"
    bam = opts.bam
    if bam and op.exists(bam):
        bam = op.abspath(bam)
        method = "GG"

    pf = inparam.split(".")[0]
    tfolder = "{0}_{1}".format(pf, method)

    cwd = os.getcwd()
    mkdir(tfolder)
    os.chdir(tfolder)

    cmds = []

    # set TRINITY_HOME env variable when preparing shell script
    env_cmd = 'export TRINITY_HOME="{0}"'.format(trinity_home)
    cmds.append(env_cmd)

    if method == "DN":
        assert op.exists("../" + inparam)

        flist = iglob("../" + inparam, opts.names)
        if paired:
            f1 = [
                x for x in flist if "_1_" in x or ".1." in x or "_1." in x or "_R1" in x
            ]
            f2 = [
                x for x in flist if "_2_" in x or ".2." in x or "_2." in x or "_R2" in x
            ]
            assert len(f1) == len(f2)
            if merge:
                r1, r2 = "left.fastq", "right.fastq"
                reads = ((f1, r1), (f2, r2))
        else:
            if merge:
                r = "single.fastq"
                reads = ((flist, r),)

        if merge:
            for fl, r in reads:
                fm = FileMerger(fl, r)
                fm.merge(checkexists=True)

    cmd = op.join(trinity_home, "Trinity")
    cmd += " --seqType fq --max_memory {0} --CPU {1}".format(opts.max_memory, opts.cpus)
    cmd += " --min_contig_length {0}".format(opts.min_contig_length)

    if opts.bflyGCThreads:
        cmd += " --bflyGCThreads {0}".format(opts.bflyGCThreads)

    if method == "GG":
        cmd += " --genome_guided_bam {0}".format(bam)
        cmd += " --genome_guided_max_intron {0}".format(opts.max_intron)
    else:
        if paired:
            if merge:
                cmd += " --left {0} --right {1}".format(reads[0][-1], reads[1][-1])
            else:
                cmd += " --left {0}".format(",".join(f1))
                cmd += " --right {0}".format(",".join(f2))
        else:
            if merge:
                cmd += " --single {0}".format(reads[0][-1])
            else:
                for f in flist:
                    cmd += " --single {0}".format(f)

    if opts.grid and opts.grid_conf_file:
        hpc_grid_runner = op.join(hpc_grid_runner_home, "hpc_cmds_GridRunner.pl")
        hpc_grid_conf_file = op.join(
            hpc_grid_runner_home, "hpc_conf", opts.grid_conf_file
        )
        assert op.exists(
            hpc_grid_conf_file
        ), "HpcGridRunner conf file does not exist: {0}".format(hpc_grid_conf_file)

        cmd += ' --grid_exec "{0} --grid_conf {1} -c"'.format(
            hpc_grid_runner, hpc_grid_conf_file
        )

    if opts.extra:
        cmd += " {0}".format(opts.extra)

    cmds.append(cmd)

    if opts.cleanup:
        cleanup_cmd = (
            'rm -rf !("Trinity.fasta"|"Trinity.gene_trans_map"|"Trinity.timing")'
            if method == "DN"
            else 'rm -rf !("Trinity-GG.fasta"|"Trinity-GG.gene_trans_map"|"Trinity.timing")'
        )
        cmds.append(cleanup_cmd)

    runfile = "run.sh"
    write_file(runfile, "\n".join(cmds))
    os.chdir(cwd)


if __name__ == "__main__":
    main()
