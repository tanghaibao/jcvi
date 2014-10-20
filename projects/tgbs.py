#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Reference-free tGBS related functions.
"""

import os.path as op
import sys

from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh


def main():

    actions = (
        ('snp', 'run SNP calling on GSNAP output'),
        ('bam', 'convert GSNAP output to BAM'),
        ('novo', 'reference-free tGBS pipeline'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def novo(args):
    """
    %prog novo reads.fastq

    Reference-free tGBS pipeline.
    """
    from jcvi.assembly.kmer import jellyfish, histogram
    from jcvi.apps.grid import MakeManager
    from jcvi.apps.cdhit import deduplicate

    p = OptionParser(novo.__doc__)
    p.add_option("--technology", choices=("illumina", "454", "iontorrent"),
                 default="iontorrent", help="Sequencing platform")
    p.set_home("cdhit")
    p.set_home("fermi")
    p.set_home("fiona")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    cpus = opts.cpus
    pf, sf = fastqfile.rsplit(".", 1)
    jf = pf + "-K23.histogram"
    if need_update(fastqfile, jf):
        jellyfish([fastqfile, "--prefix={0}".format(pf),
                    "--cpus={0}".format(cpus)])

    genomesize = histogram([jf, pf, "23"])
    fiona = pf + ".fiona." + sf
    if need_update(fastqfile, fiona):
        cmd = op.join(opts.fiona_home, "bin/fiona")
        cmd += " -g {0} -nt {1} --sequencing-technology {2}".\
                    format(genomesize, cpus, opts.technology)
        cmd += " -vv {0} {1}".format(fastqfile, fiona)
        logfile = pf + ".fiona.log"
        sh(cmd, outfile=logfile, errfile=logfile)

    mag = pf + ".p0.mag.gz"
    fhome = opts.fermi_home
    if need_update(fiona, mag):
        cmd = op.join(fhome, "run-fermi.pl")
        cmd += " -t {0} -e {1}/fermi -p {2} {3}".format(cpus, fhome, pf, fiona)
        sh(cmd, outfile="makefile")
        MakeManager().run(cpus=cpus)

    asm = pf + ".fermi.fasta"
    if need_update(mag, asm):
        cmd = "seqtk seq {0} -A".format(mag)
        sh(cmd, outfile=asm)

    cons = asm + ".P98.cdhit.consensus.fasta"
    if need_update(asm, cons):
        deduplicate([asm, "--consensus", "--reads",
                      "--cdhit_home={0}".format(opts.cdhit_home)])


def bam(args):
    """
    %prog snp input.gsnap ref.fasta

    Convert GSNAP output to BAM.
    """
    from jcvi.formats.sizes import Sizes
    from jcvi.formats.sam import index

    p = OptionParser(bam.__doc__)
    p.set_home("eddyyeh")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gsnapfile, fastafile = args
    EYHOME = opts.eddyyeh_home
    pf = gsnapfile.rsplit(".", 1)[0]
    uniqsam = pf + ".unique.sam"
    if need_update((gsnapfile, fastafile), uniqsam):
        cmd = op.join(EYHOME, "gsnap2gff3.pl")
        sizesfile = Sizes(fastafile).filename
        cmd += " --format sam -i {0} -o {1}".format(gsnapfile, uniqsam)
        cmd += " -u -l {0} -p {1}".format(sizesfile, opts.cpus)
        sh(cmd)

    index([uniqsam])


def snp(args):
    """
    %prog snp input.gsnap

    Run SNP calling on GSNAP output after apps.gsnap.align().
    """
    p = OptionParser(snp.__doc__)
    p.set_home("eddyyeh")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gsnapfile, = args
    EYHOME = opts.eddyyeh_home
    pf = gsnapfile.rsplit(".", 1)[0]
    nativefile = pf + ".native"
    if need_update(gsnapfile, nativefile):
        cmd = op.join(EYHOME, "convert2native.pl")
        cmd += " --gsnap {0} -o {1}".format(gsnapfile, nativefile)
        cmd += " -proc {0}".format(opts.cpus)
        sh(cmd)

    snpfile = pf + ".snp"
    if need_update(nativefile, snpfile):
        cmd = op.join(EYHOME, "SNPs/SNP_Discovery-short.pl")
        cmd += " --native {0} -o {1}".format(nativefile, snpfile)
        cmd += " -a 2 -ac 0.3 -c 0.8"
        sh(cmd)


if __name__ == '__main__':
    main()
