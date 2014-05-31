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
    p.set_home("trinity")
    p.set_cpus()
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

    cmd = op.join(thome, "Trinity.pl")
    cmd += " --seqType fq --JM 100G --CPU {0}".format(opts.cpus)
    if paired:
        cmd += " --left {0} --right {1}".format(reads[0][-1], reads[1][-1])
    else:
        cmd += " --single {0}".format(reads[0][-1])

    runfile = "run.sh"
    write_file(runfile, cmd, meta="run script")
    os.chdir(cwd)


if __name__ == '__main__':
    main()
