#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper to trim and correct sequence data.
"""

import os
import os.path as op
import sys

from optparse import OptionParser

from jcvi.formats.fastq import guessoffset, convert
from jcvi.apps.command import JAVA
from jcvi.apps.base import ActionDispatcher, debug, set_grid, download, sh
debug()


def main():

    actions = (
        ('trim', 'trim reads using TRIMMOMATIC'),
        ('correct', 'correct reads with ALLPATHS-LG'),
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
    p = OptionParser(trim.__doc__)
    p.add_option("--path", default="/home/htang/bin/trimmomatic-0.13.jar",
            help="Path to trimmomatic [default: %default]")
    set_grid(p)

    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    path = opts.path
    url = \
    "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.13.zip"

    if not op.exists(path):
        path = download(url)
        if not op.exists("Trimmomatic-0.13"):
            sh("unzip " + path)
        os.remove(path)
        path = "Trimmomatic-0.13/trimmomatic-0.13.jar"

    assert op.exists(path)

    adaptersfile = "adapters.fasta"
    if not op.exists(adaptersfile):
        sh("cp ~/adapters.fasta .")

    assert op.exists("adapters.fasta"), \
        "Please place the illumina adapter sequence in `adapters.fasta`"

    # Trimmomatic only accepts offset-64, convert if needed
    for i, fastqfile in enumerate(args):
        offset = guessoffset([fastqfile])
        if offset != 64:
            newfastqfile = fastqfile.rstrip(".gz").rsplit(".", 1)[0] + ".b64.fastq.gz"
            convert([fastqfile, newfastqfile])
            args[i] = newfastqfile

    cmd = op.join(JAVA, "java-1.6.0 -Xmx4g")
    cmd += " -cp {0} org.usadellab.trimmomatic".format(path)
    frags = ".frags.fastq.gz"
    pairs = ".pairs.fastq.gz"
    if len(args) == 1:
        cmd += ".TrimmomaticSE"
        fastqfile, = args
        prefix = fastqfile.split(".")[0]
        frags1 = prefix + frags
        cmd += " {0}".format(" ".join((fastqfile, frags1)))
    else:
        cmd += ".TrimmomaticPE"
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
    cmd += " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    sh(cmd, grid=opts.grid)


if __name__ == '__main__':
    main()
