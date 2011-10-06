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
from jcvi.apps.base import ActionDispatcher, debug, set_grid, download, sh
debug()


def main():

    actions = (
        ('trim', 'trim reads using TRIMMOMATIC'),
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


if __name__ == '__main__':
    main()
