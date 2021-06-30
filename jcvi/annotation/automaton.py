#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Automate genome annotation by iterating processing a set of files, individually.
"""

import os.path as op
import shutil
import sys
import logging

from functools import partial
from tempfile import mkdtemp

from jcvi.assembly.automaton import iter_project
from jcvi.apps.grid import Jobs, MakeManager
from jcvi.formats.base import FileMerger, split
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, mkdir, sh, iglob


def main():

    actions = (
        ("augustus", "run parallel AUGUSTUS"),
        ("cufflinks", "run cufflinks following tophat"),
        ("star", "run star alignment"),
        ("tophat", "run tophat on a list of inputs"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def augustuswrap(fastafile, species="maize", gff3=True, cfgfile=None, hintsfile=None):
    cmd = "augustus {0}".format(fastafile)
    if gff3:
        cmd += " --gff3=on"
    cmd += " --species={0}".format(species)
    if cfgfile:
        cmd += " --extrinsicCfgFile={0}".format(cfgfile)
    if hintsfile:
        cmd += " --alternatives-from-evidence=true"
        cmd += " --hintsfile={0} --allow_hinted_splicesites=atac".format(hintsfile)
    cmd += " --introns=on --genemodel=complete"
    suffix = ".gff3" if gff3 else ".out"
    outfile = fastafile.rsplit(".", 1)[0] + suffix
    sh(cmd, outfile=outfile)
    return outfile


def augustus(args):
    """
    %prog augustus fastafile

    Run parallel AUGUSTUS. Final results can be reformatted using
    annotation.reformat.augustus().
    """
    p = OptionParser(augustus.__doc__)
    p.add_option("--species", default="maize", help="Use species model for prediction")
    p.add_option("--hintsfile", help="Hint-guided AUGUSTUS")
    p.add_option("--nogff3", default=False, action="store_true", help="Turn --gff3=off")
    p.set_home("augustus")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    cpus = opts.cpus
    mhome = opts.augustus_home
    gff3 = not opts.nogff3
    suffix = ".gff3" if gff3 else ".out"
    cfgfile = op.join(mhome, "config/extrinsic/extrinsic.M.RM.E.W.cfg")

    outdir = mkdtemp(dir=".")
    fs = split([fastafile, outdir, str(cpus)])

    augustuswrap_params = partial(
        augustuswrap,
        species=opts.species,
        gff3=gff3,
        cfgfile=cfgfile,
        hintsfile=opts.hintsfile,
    )
    g = Jobs(augustuswrap_params, fs.names)
    g.run()

    gff3files = [x.rsplit(".", 1)[0] + suffix for x in fs.names]
    outfile = fastafile.rsplit(".", 1)[0] + suffix
    FileMerger(gff3files, outfile=outfile).merge()
    shutil.rmtree(outdir)

    if gff3:
        from jcvi.annotation.reformat import augustus as reformat_augustus

        reformat_outfile = outfile.replace(".gff3", ".reformat.gff3")
        reformat_augustus([outfile, "--outfile={0}".format(reformat_outfile)])


def star(args):
    """
    %prog star folder reference

    Run star on a folder with reads.
    """
    p = OptionParser(star.__doc__)
    p.add_option(
        "--single", default=False, action="store_true", help="Single end mapping"
    )
    p.set_fastq_names()
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

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
    p.add_option("--gtf", help="Reference annotation")
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

    mergedgtf = "merged/merged.gtf"
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
    from jcvi.apps.bowtie import check_index
    from jcvi.formats.fastq import guessoffset

    p = OptionParser(tophat.__doc__)
    p.add_option("--gtf", help="Reference annotation")
    p.add_option(
        "--single", default=False, action="store_true", help="Single end mapping"
    )
    p.add_option(
        "--intron",
        default=15000,
        type="int",
        help="Max intron size",
    )
    p.add_option(
        "--dist",
        default=-50,
        type="int",
        help="Mate inner distance",
    )
    p.add_option(
        "--stdev",
        default=50,
        type="int",
        help="Mate standard deviation",
    )
    p.set_phred()
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    num = 1 if opts.single else 2
    folder, reference = args
    reference = check_index(reference)
    for p, prefix in iter_project(folder, n=num):
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
            (a,) = p
        else:  # Paired-end
            a, b = p
            cmd += " --max-intron-length {0}".format(opts.intron)
            cmd += " --mate-inner-dist {0}".format(opts.dist)
            cmd += " --mate-std-dev {0}".format(opts.stdev)

        phred = opts.phred or str(guessoffset([a]))
        if phred == "64":
            cmd += " --phred64-quals"
        cmd += " {0} {1}".format(reference, " ".join(p))

        sh(cmd)


if __name__ == "__main__":
    main()
