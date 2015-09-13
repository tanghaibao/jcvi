#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Automate genome annotation by iterating processing a set of files, individually.
"""

import os
import os.path as op
import sys
import logging

from jcvi.assembly.automaton import iter_project
from jcvi.apps.grid import MakeManager
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, \
            mkdir, sh, iglob


def main():

    actions = (
        ('cufflinks', 'run cufflinks following tophat'),
        ('star', 'run star alignment'),
        ('tophat', 'run tophat on a list of inputs'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def star(args):
    """
    %prog star folder reference

    Run star on a folder with reads.
    """
    p = OptionParser(star.__doc__)
    p.add_option("--single", default=False, action="store_true",
                 help="Single end mapping")
    p.set_fastq_names()
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, reference = args
    cpus = opts.cpus
    mm = MakeManager()

    num = 1 if opts.single else 2
    folder, reference = args
    gd = "GenomeDir"
    mkdir(gd)
    STAR = "STAR --runThreadN {0} --genomeDir {1}".format(cpus, gd)

    # Step 0: build genome index
    genomeidx = op.join(gd, "Genome")
    if need_update(reference, genomeidx):
        cmd = STAR + " --runMode genomeGenerate"
        cmd += " --genomeFastaFiles {0}".format(reference)
        mm.add(reference, genomeidx, cmd)

    # Step 1: align
    for p, prefix in iter_project(folder, opts.names, num):
        pf = "{0}_star".format(prefix)
        bamfile = pf + "Aligned.sortedByCoord.out.bam"
        cmd = STAR + " --readFilesIn {0}".format(" ".join(p))
        if p[0].endswith(".gz"):
            cmd += " --readFilesCommand zcat"
        cmd += " --outSAMtype BAM SortedByCoordinate"
        cmd += " --outFileNamePrefix {0}".format(pf)
        cmd += " --twopassMode Basic"
        # Compatibility for cufflinks
        cmd += " --outSAMstrandField intronMotif"
        cmd += " --outFilterIntronMotifs RemoveNoncanonical"
        mm.add(p, bamfile, cmd)

    mm.write()


def cufflinks(args):
    """
    %prog cufflinks folder reference

    Run cufflinks on a folder containing tophat results.
    """
    p = OptionParser(cufflinks.__doc__)
    p.add_option("--gtf", help="Reference annotation [default: %default]")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, reference = args
    cpus = opts.cpus
    gtf = opts.gtf
    transcripts = "transcripts.gtf"

    mm = MakeManager()
    gtfs = []
    for bam in iglob(folder, "*.bam"):
        pf = op.basename(bam).split(".")[0]
        outdir = pf + "_cufflinks"
        cmd = "cufflinks"
        cmd += " -o {0}".format(outdir)
        cmd += " -p {0}".format(cpus)
        if gtf:
            cmd += " -g {0}".format(gtf)
        cmd += " --frag-bias-correct {0}".format(reference)
        cmd += " --multi-read-correct"
        cmd += " {0}".format(bam)
        cgtf = op.join(outdir, transcripts)
        mm.add(bam, cgtf, cmd)
        gtfs.append(cgtf)

    assemblylist = "assembly_list.txt"
    cmd = 'find . -name "{0}" > {1}'.format(transcripts, assemblylist)
    mm.add(gtfs, assemblylist, cmd)

    mergedgtf = "merged.gtf"
    cmd = "cuffmerge"
    cmd += " -o merged"
    cmd += " -p {0}".format(cpus)
    if gtf:
        cmd += " -g {0}".format(gtf)
    cmd += " -s {0}".format(reference)
    cmd += " {0}".format(assemblylist)
    mm.add(assemblylist, mergedgtf, cmd)

    mm.write()


def tophat(args):
    """
    %prog tophat folder reference

    Run tophat on a folder of reads.
    """
    p = OptionParser(tophat.__doc__)
    p.add_option("--gtf", help="Reference annotation [default: %default]")
    p.add_option("--single", default=False, action="store_true",
                 help="Single end mapping")
    p.add_option("--intron", default=15000, type="int",
                 help="Max intron size [default: %default]")
    p.add_option("--dist", default=-50, type="int",
                 help="Mate inner distance [default: %default]")
    p.add_option("--stdev", default=50, type="int",
                 help="Mate standard deviation [default: %default]")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    num = 1 if opts.single else 2
    folder, reference = args
    for p, prefix in iter_project(folder, num):
        outdir = "{0}_tophat".format(prefix)
        outfile = op.join(outdir, "accepted_hits.bam")
        if op.exists(outfile):
            logging.debug("File `{0}` found. Skipping.".format(outfile))
            continue

        cmd = "tophat -p {0}".format(opts.cpus)
        if opts.gtf:
            cmd += " -G {0}".format(opts.gtf)
        cmd += " -o {0}".format(outdir)

        if num == 1:  # Single-end
            a, = p
            cmd += " {0} {1}".format(reference, a)
        else:  # Paired-end
            a, b = p
            cmd += " --max-intron-length {0}".format(opts.intron)
            cmd += " --mate-inner-dist {0}".format(opts.dist)
            cmd += " --mate-std-dev {0}".format(opts.stdev)
            cmd += " {0} {1} {2}".format(reference, a, b)

        sh(cmd)


if __name__ == '__main__':
    main()
