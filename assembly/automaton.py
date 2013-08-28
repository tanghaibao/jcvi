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
from jcvi.apps.base import MOptionParser

from jcvi.utils.iter import grouper
from jcvi.apps.base import ActionDispatcher, debug, need_update, mkdir, sh
debug()


def main():

    actions = (
        ('correct', 'run automated ALLPATHS correction'),
        ('allpaths', 'run automated ALLPATHS'),
        ('soap', 'run automated SOAP'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def slink(p, pf, tag, extra=None):

    mkdir(pf, overwrite=True)
    cwd = os.getcwd()
    os.chdir(pf)

    # Create sym-links for the input files
    i = 1
    for f in sorted(p):
        gz = ".gz" if f.endswith(".gz") else ""
        if "PE-0" in f:
            sh("ln -sf ../{0} PE-0.fastq{1}".format(f, gz))
            continue
        for t in tag:
            sh("ln -sf ../{0} {1}.{2}.fastq{3}".format(f, t, i, gz))
        i += 1

    if extra:
        for e in extra:
            sh("ln -sf {0}".format(e))

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
    prepare([pf] + sorted(glob("*.fastq") + glob("*.fastq.gz")))
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
    targets = (cm + ".1.corr.fastq", cm + ".2.corr.fastq", \
                pf + ".PE-0.corr.fastq")
    if not need_update(p, targets):
        logging.debug("Corrected reads found: {0}. Skipped.".format(targets))
        return

    slink(p, pf, tag)

    cwd = os.getcwd()
    os.chdir(pf)
    cr(sorted(glob("*.fastq") + glob("*.fastq.gz")) + ["--nofragsdedup"])
    sh("mv {0}.1.corr.fastq ../{1}".format(itag, targets[0]))
    sh("mv {0}.2.corr.fastq ../{1}".format(itag, targets[1]))
    sh("mv frag_reads_corr.corr.fastq ../{0}".format(targets[2]))

    logging.debug("Correction finished: {0}".format(targets))
    os.chdir(cwd)


def soap_trios(p, pf, tag, extra):
    """
    Take one pair of reads and 'widow' reads after correction and run SOAP.
    """
    from jcvi.assembly.soap import prepare

    logging.debug("Work on {0} ({1})".format(pf, ','.join(p)))
    asm = "{0}.closed.scafSeq".format(pf)
    if not need_update(p, asm):
        logging.debug("Assembly found: {0}. Skipped.".format(asm))
        return

    slink(p, pf, tag, extra)

    cwd = os.getcwd()
    os.chdir(pf)
    prepare(sorted(glob("*.fastq") + glob("*.fastq.gz")) + \
            ["--assemble_1st_rank_only", "-K 31"])
    sh("./run.sh")
    sh("cp asm31.closed.scafSeq ../{0}".format(asm))

    logging.debug("Assembly finished: {0}".format(asm))
    os.chdir(cwd)


def iter_project(folder, n=2):
    # Check for paired reads and extract project id
    filelist = sorted(glob(folder + "/*.*"))
    for p in grouper(n, filelist):
        if len(p) != n:
            continue

        pp = [op.basename(x) for x in p]
        pf = op.commonprefix(pp).strip("._-")
        yield p, pf


def soap(args):
    """
    %prog soap folder tag [*.fastq]

    Run SOAP on a folder of paired reads and apply tag before assembly.
    Optional *.fastq in the argument list will be symlinked in each folder and
    co-assembled.
    """
    from jcvi.apps.softlink import get_abs_path

    p = MOptionParser(soap.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    folder, tag = args[:2]
    extra = args[2:]
    extra = [get_abs_path(x) for x in extra]
    tag = tag.split(",")
    for p, pf in iter_project(folder, n=3):
        soap_trios(p, pf, tag, extra)


def correct(args):
    """
    %prog correct folder tag

    Run ALLPATHS correction on a folder of paired reads and apply tag.
    """
    p = MOptionParser(correct.__doc__)
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
    p = MOptionParser(allpaths.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, tag = args
    tag = tag.split(",")
    for p, pf in iter_project(folder):
        assemble_pairs(p, pf, tag)


if __name__ == '__main__':
    main()
