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
from jcvi.apps.command import JAVAPATH
from jcvi.apps.base import ActionDispatcher, debug, set_grid, download, \
        sh, mkdir, need_update
debug()


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
    TrimVersion = tv = "0.15"
    TrimJar = "trimmomatic-{0}.jar".format(tv)
    p = OptionParser(trim.__doc__)
    p.add_option("--path", default=op.join("~/bin", TrimJar),
            help="Path to trimmomatic [default: %default]")
    p.add_option("--phred", default=None, type="int",
            help="Phred score offset [default: %default]")
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
        sh("cp ~/adapters.fasta .")

    assert op.exists("adapters.fasta"), \
        "Please place the illumina adapter sequence in `adapters.fasta`"

    if opts.phred is None:
        offset = guessoffset([args[0]])
    else:
        offset = opts.phred

    phredflag = " -phred{0}".format(offset)

    cmd = JAVAPATH("java-1.6.0")
    cmd += " -Xmx4g -cp {0} org.usadellab.trimmomatic".format(path)
    frags = ".frags.fastq.gz"
    pairs = ".pairs.fastq.gz"
    if len(args) == 1:
        cmd += ".TrimmomaticSE"
        cmd += phredflag
        fastqfile, = args
        prefix = fastqfile.replace(".gz", "").rsplit(".", 1)[0]
        frags1 = prefix + frags
        cmd += " {0}".format(" ".join((fastqfile, frags1)))
    else:
        cmd += ".TrimmomaticPE"
        cmd += phredflag
        fastqfile1, fastqfile2 = args
        prefix1 = fastqfile1.split(".")[0]
        prefix2 = fastqfile2.split(".")[0]
        pairs1 = prefix1 + pairs
        pairs2 = prefix2 + pairs
        frags1 = prefix1 + frags
        frags2 = prefix2 + frags
        cmd += " {0}".format(" ".join((fastqfile1, fastqfile2, \
                pairs1, frags1, pairs2, frags2)))

    cmd += " ILLUMINACLIP:adapters.fasta:2:40:15"
    cmd += " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"
    sh(cmd, grid=opts.grid)


def correct(args):
    """
    %prog correct fastqfile

    Correct the fastqfile and generated corrected fastqfiles. If --paired is
    used, then will write to `pairs.corr.fastq` and `frags.corr.fastq`. Remember
    the pairs need to be in interleaved format.

    Rename your fastqfile to look like the following:
    PE-376.fastq (paired end)
    MP-3000.fastq (mate pairs)
    TT-3000.fastq (mate pairs, but from 454 data)

    The insert size does NOT affect correction, but MP- reads will be removed of
    duplicates.
    """
    from jcvi.assembly.allpaths import prepare

    p = OptionParser(correct.__doc__)
    p.add_option("--paired", default=False, action="store_true",
                 help="Input is interleaved fastq file [default: %default]")
    p.add_option("--cpus", default=32, type="int",
                 help="Number of threads to run [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastq, = args
    prefix = op.basename(fastq).rsplit(".", 1)[0]
    mplib = fastq.startswith("MP-")
    tag = "jump" if mplib else "frag"
    tag += "_reads"

    prepare(["Unknown", fastq])
    sh("ulimit -s 100000")

    datadir = "data"
    mkdir(datadir)
    fullpath = op.join(os.getcwd(), datadir)
    nthreads = " NUM_THREADS={0}".format(opts.cpus)

    assert op.exists(fastq)
    orig = datadir + "/{0}_orig".format(tag)
    origfastb = orig + ".fastb"
    if need_update(fastq, origfastb):
        cmd = "PrepareAllPathsInputs.pl DATA_DIR={0}".format(fullpath)
        sh(cmd)

    assert op.exists(origfastb)
    filt = datadir + "/{0}_filt".format(tag)
    filtfastb = filt + ".fastb"
    if need_update(origfastb, filtfastb):
        cmd = "RemoveDodgyReads IN_HEAD={0} OUT_HEAD={1}".format(orig, filt)
        if fastq.startswith("MP-"):
            cmd += " REMOVE_DUPLICATES=True"
        cmd += nthreads
        sh(cmd)

    assert op.exists(filtfastb)
    prec = datadir + "/{0}_prec".format(tag)
    precfastb = prec + ".fastb"
    if need_update(filtfastb, precfastb):
        cmd = "PreCorrect IN_HEAD={0} OUT_HEAD={1}".format(filt, prec)
        cmd += nthreads
        sh(cmd)

    filtpairs = filt + ".pairs"
    precpairs = prec + ".pairs"
    if need_update(filtpairs, precpairs):
        cmd = "ln -sf {0}.pairs {1}.pairs".format(op.basename(filt), prec)
        sh(cmd)

    assert op.exists(precfastb)
    edit = datadir + "/{0}_edit".format(tag)
    editfastb = edit + ".fastb"
    corr = datadir + "/{0}_corr".format(tag)
    corrfastb = corr + ".fastb"
    if need_update(precfastb, corrfastb):
        cmd = "FindErrors DELETE=True IN_HEAD={0}".format(prec)
        cmd += " OUT_EDIT_HEAD={0} OUT_CORR_HEAD={1}".format(edit, corr)
        sh(cmd)

    assert op.exists(corrfastb)
    corrfastq = corr + ".fastq"
    if need_update(corrfastb, corrfastq):
        cmd = "FastbAndQualb2Fastq HEAD={0}".format(corr)
        sh(cmd)

    assert op.exists(corrfastq)


if __name__ == '__main__':
    main()
