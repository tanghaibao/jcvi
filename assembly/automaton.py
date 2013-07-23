#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Automate genome assembly by iterating assembly on a set of files, individually.
This is useful for BAC pool assembly where you want to assemble per pool.
"""

import os
import os.path as op
import sys
import logging

from glob import glob
from optparse import OptionParser

from jcvi.utils.iter import grouper
from jcvi.apps.base import ActionDispatcher, debug, need_update, mkdir, sh
debug()


def main():

    actions = (
        ('correct', 'run automated ALLPATHS correction'),
        ('allpaths', 'run automated ALLPATHS'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def slink(p, pf, tag):
    mkdir(pf, overwrite=True)
    cwd = os.getcwd()
    os.chdir(pf)
    gz = ".gz" if p[0].endswith(".gz") else ""

    # Create sym-links for the input files
    for i, f in enumerate(sorted(p)):
        for t in tag:
            sh("ln -sf ../{0} {1}.{2}.fastq{3}".format(f, t, i + 1, gz))
    os.chdir(cwd)


def assemble_pairs(p, pf, tag):
    """
    Take one pair of reads and assemble to contigs.fasta.
    """
    from jcvi.assembly.preprocess import prepare

    logging.debug("Work on {0} ({1})".format(pf, ','.join(p)))
    asm = "{0}.contigs.fasta".format(pf)
    if not need_update(p, asm):
        logging.debug("Assembly found: {0}. Skipped.".format(asm))
        return

    slink(p, pf, tag)

    cwd = os.getcwd()
    os.chdir(pf)
    prepare([pf] + glob("*.fastq") + glob("*.fastq.gz"))
    sh("chmod u+x run.sh")
    sh("./run.sh")
    sh("cp allpaths/ASSEMBLIES/run/final.contigs.fasta ../{0}".format(asm))

    logging.debug("Assembly finished: {0}".format(asm))
    os.chdir(cwd)


def correct_pairs(p, pf, tag):
    """
    Take one pair of reads and correct to generate *.corr.fastq.
    """
    from jcvi.assembly.preprocess import correct as cr

    logging.debug("Work on {0} ({1})".format(pf, ','.join(p)))
    itag = tag[0]
    cm = ".".join((pf, itag))
    targets = (cm + ".1.corr.fastq", cm + ".2.corr.fastq")
    if not need_update(p, targets):
        logging.debug("Corrected reads found: {0}. Skipped.".format(targets))
        return

    slink(p, pf, tag)

    cwd = os.getcwd()
    os.chdir(pf)
    cr(glob("*.fastq") + glob("*.fastq.gz"))
    sh("mv {0}.1.corr.fastq ../{1}".format(itag, targets[0]))
    sh("mv {0}.2.corr.fastq ../{1}".format(itag, targets[1]))

    logging.debug("Correction finished: {0}".format(targets))
    os.chdir(cwd)


def iter_project(folder):
    # Check for paired reads and extract project id
    filelist = sorted(glob(folder + "/*.*"))
    for p in grouper(2, filelist):
        if len(p) != 2:
            continue

        pp = [op.basename(x) for x in p]
        pf = op.commonprefix(pp).strip("._-")
        yield p, pf


def correct(args):
    """
    %prog correct folder tag

    Run ALLPATHS correction on a folder of paired reads and apply tag.
    """
    p = OptionParser(correct.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, tag = args
    tag = tag.split(",")
    for p, pf in iter_project(folder):
        correct_pairs(p, pf, tag)


def allpaths(args):
    """
    %prog automaton folder tag

    Run assembly on a folder of paired reads and apply tag (PE-200, PE-500).
    Allow multiple tags separated by comma, e.g. PE-350,TT-1050
    """
    p = OptionParser(allpaths.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, tag = args
    tag = tag.split(",")
    for p, pf in iter_project(folder):
        assemble_pairs(p, pf, tag)


if __name__ == '__main__':
    main()
