#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper to trim and correct sequence data.
"""

import os
import os.path as op
import sys

from optparse import OptionParser

from jcvi.formats.fastq import guessoffset
from jcvi.utils.cbook import depends
from jcvi.apps.command import JAVAPATH
from jcvi.apps.base import ActionDispatcher, debug, set_grid, download, \
        sh, mkdir, write_file, need_update
debug()


Adapters = """
>Illumina_PE-1rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>Illumina_PE-2rc
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
"""


def main():

    actions = (
        ('trim', 'trim reads using TRIMMOMATIC'),
        ('correct', 'correct reads using ALLPATHS-LG'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def trim(args):
    """
    %prog trim fastqfiles

    Trim reads using TRIMMOMATIC. If two fastqfiles are given, then it invokes
    the paired reads mode. See manual:

    <http://www.usadellab.org/cms/index.php?page=trimmomatic>
    """
    TrimVersion = tv = "0.20"
    TrimJar = "trimmomatic-{0}.jar".format(tv)
    phdchoices = ("33", "64")
    p = OptionParser(trim.__doc__)
    p.add_option("--path", default=op.join("~/bin", TrimJar),
            help="Path to trimmomatic [default: %default]")
    p.add_option("--phred", default=None, choices=phdchoices,
            help="Phred score offset {0} [default: guess]".format(phdchoices))
    p.add_option("--nofrags", default=False, action="store_true",
            help="Discard frags file in PE mode [default: %default]")
    p.add_option("--minqv", default=10, type="int",
            help="Average qv after trimming [default: %default]")
    p.add_option("--minlen", default=30, type="int",
            help="Minimum length after trimming [default: %default]")
    p.add_option("--nogz", default=False, action="store_true",
            help="Do not write to gzipped files [default: %default]")
    set_grid(p)

    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    path = op.expanduser(opts.path)
    url = \
    "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-{0}.zip"\
    .format(tv)

    if not op.exists(path):
        path = download(url)
        TrimUnzipped = "Trimmomatic-" + tv
        if not op.exists(TrimUnzipped):
            sh("unzip " + path)
        os.remove(path)
        path = op.join(TrimUnzipped, TrimJar)

    assert op.exists(path)

    adaptersfile = "adapters.fasta"
    if not op.exists(adaptersfile):
        write_file(adaptersfile, Adapters)

    assert op.exists(adaptersfile), \
        "Please place the illumina adapter sequence in `{0}`".\
        format(adaptersfile)

    if opts.phred is None:
        offset = guessoffset([args[0]])
    else:
        offset = int(opts.phred)

    phredflag = " -phred{0}".format(offset)

    cmd = JAVAPATH("java-1.6.0")
    cmd += " -Xmx4g -cp {0} org.usadellab.trimmomatic".format(path)
    frags = ".frags.fastq"
    pairs = ".pairs.fastq"
    if not opts.nogz:
        frags += ".gz"
        pairs += ".gz"

    get_prefix = lambda x: op.basename(x).replace(".gz", "").rsplit(".", 1)[0]
    if len(args) == 1:
        cmd += ".TrimmomaticSE"
        cmd += phredflag
        fastqfile, = args
        prefix = get_prefix(fastqfile)
        frags1 = prefix + frags
        cmd += " {0}".format(" ".join((fastqfile, frags1)))
    else:
        cmd += ".TrimmomaticPE"
        cmd += phredflag
        fastqfile1, fastqfile2 = args
        prefix1 = get_prefix(fastqfile1)
        prefix2 = get_prefix(fastqfile2)
        pairs1 = prefix1 + pairs
        pairs2 = prefix2 + pairs
        frags1 = prefix1 + frags
        frags2 = prefix2 + frags
        if opts.nofrags:
            frags1 = "/dev/null"
            frags2 = "/dev/null"
        cmd += " {0}".format(" ".join((fastqfile1, fastqfile2, \
                pairs1, frags1, pairs2, frags2)))

    cmd += " ILLUMINACLIP:{0}:2:40:12".format(adaptersfile)
    cmd += " LEADING:3 TRAILING:3"
    cmd += " SLIDINGWINDOW:4:{0} MINLEN:{1}".format(opts.minqv, opts.minlen)
    if offset != 33:
        cmd += " TOPHRED33"
    sh(cmd, grid=opts.grid)


@depends
def run_RemoveDodgyReads(infile=None, outfile=None, workdir=None,
        removeDuplicates=True, rc=False, nthreads=32):
    # orig.fastb => filt.fastb
    assert op.exists(infile)
    orig = infile.rsplit(".", 1)[0]
    filt = outfile.rsplit(".", 1)[0]

    cmd = "RemoveDodgyReads IN_HEAD={0} OUT_HEAD={1}".format(orig, filt)
    if not removeDuplicates:
        cmd += " REMOVE_DUPLICATES=False"
    if rc:
        cmd += " RC=True"
    cmd += nthreads
    sh(cmd)


@depends
def run_FastbAndQualb2Fastq(infile=None, outfile=None):
    corr = op.basename(infile).rsplit(".", 1)[0]
    cmd = "FastbQualbToFastq HEAD_IN={0} HEAD_OUT={0}".format(corr)
    cmd += " PAIRED=False PHRED_OFFSET=33"
    sh(cmd)


@depends
def run_pairs(infile=None, outfile=None):
    from jcvi.assembly.allpaths import pairs
    pairs(infile)


def correct(args):
    """
    %prog correct *.fastq

    Correct the fastqfile and generated corrected fastqfiles. This calls
    assembly.allpaths.prepare() to generate input files for ALLPATHS-LG. The
    naming convention for your fastqfiles are important, and are listed below.

    By default, this will correct all PE reads, and remove duplicates of all MP
    reads, and results will be placed in `frag_reads.corr.{pairs,frags}.fastq`
    and `jump_reads.corr.{pairs,frags}.fastq`.
    """
    from jcvi.assembly.allpaths import prepare
    from jcvi.assembly.base import FastqNamings

    p = OptionParser(correct.__doc__ + FastqNamings)
    p.add_option("--nofragsdedup", default=False, action="store_true",
                 help="Don't deduplicate the fragment reads [default: %default]")
    p.add_option("--cpus", default=32, type="int",
                 help="Number of threads to run [default: %default]")
    p.add_option("--phred64", default=False, action="store_true",
                 help="Reads are all phred 64 offset [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fastq = args
    tag, tagj = "frag_reads", "jump_reads"

    prepare(["Unknown"] + fastq + ["--norun"])

    datadir = "data"
    mkdir(datadir)
    fullpath = op.join(os.getcwd(), datadir)
    nthreads = " NUM_THREADS={0}".format(opts.cpus)

    orig = datadir + "/{0}_orig".format(tag)
    origfastb = orig + ".fastb"
    if need_update(fastq, origfastb):
        cmd = "PrepareAllPathsInputs.pl DATA_DIR={0} HOSTS='{1}'".\
                format(fullpath, opts.cpus)
        if opts.phred64:
            cmd += " PHRED_64=True"
        sh(cmd)

    if op.exists(origfastb):
        dedup = not opts.nofragsdedup
        correct_frag(datadir, tag, origfastb, nthreads, dedup=dedup)

    origj = datadir + "/{0}_orig".format(tagj)
    origjfastb = origj + ".fastb"

    if op.exists(origjfastb):
        correct_jump(datadir, tagj, origjfastb, nthreads)


def correct_frag(datadir, tag, origfastb, nthreads, dedup=False):
    filt = datadir + "/{0}_filt".format(tag)
    filtfastb = filt + ".fastb"
    run_RemoveDodgyReads(infile=origfastb, outfile=filtfastb,
                         removeDuplicates=dedup, rc=False, nthreads=nthreads)

    filtpairs = filt + ".pairs"
    edit = datadir + "/{0}_edit".format(tag)
    editpairs = edit + ".pairs"
    if need_update(filtpairs, editpairs):
        cmd = "ln -sf {0} {1}.pairs".format(op.basename(filtpairs), edit)
        sh(cmd)

    editfastb = edit + ".fastb"
    if need_update(filtfastb, editfastb):
        cmd = "FindErrors HEAD_IN={0} HEAD_OUT={1}".format(filt, edit)
        cmd += nthreads
        sh(cmd)

    corr = datadir + "/{0}_corr".format(tag)
    corrfastb = corr + ".fastb"
    if need_update(editfastb, corrfastb):
        cmd = "CleanCorrectedReads DELETE=True"
        cmd += " HEAD_IN={0} HEAD_OUT={1}".format(edit, corr)
        cmd += nthreads
        sh(cmd)

    pf = op.basename(corr)

    cwd = os.getcwd()
    os.chdir(datadir)
    corrfastq = pf + ".fastq"
    run_FastbAndQualb2Fastq(infile=op.basename(corrfastb), outfile=corrfastq)
    os.chdir(cwd)

    pairsfile = pf + ".pairs"
    fragsfastq = pf + ".corr.fastq"
    run_pairs(infile=[op.join(datadir, pairsfile), op.join(datadir, corrfastq)],
                      outfile=fragsfastq)


def correct_jump(datadir, tagj, origjfastb, nthreads):
    # Pipeline for jump reads does not involve correction
    filt = datadir + "/{0}_filt".format(tagj)
    filtfastb = filt + ".fastb"
    run_RemoveDodgyReads(infile=origjfastb, outfile=filtfastb, \
                         removeDuplicates=True, rc=True, nthreads=nthreads)

    pf = op.basename(filt)

    cwd = os.getcwd()
    os.chdir(datadir)
    filtfastq = pf + ".fastq"
    run_FastbAndQualb2Fastq(infile=op.basename(filtfastb), outfile=filtfastq)
    os.chdir(cwd)

    pairsfile = pf + ".pairs"
    fragsfastq = pf + ".corr.fastq"
    run_pairs(infile=[op.join(datadir, pairsfile), op.join(datadir, filtfastq)],
                      outfile=fragsfastq)


if __name__ == '__main__':
    main()
