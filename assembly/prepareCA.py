#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Prepare input files for Celera Assembler, dispatch based on file suffix::

*.fasta: convert-fasta-to-v2.pl
*.sff: sffToCA
*.fastq: fastqToCA
"""

import os.path as op
import sys
import logging

from glob import glob
from optparse import OptionParser
from collections import defaultdict

from Bio import SeqIO

from jcvi.formats.fasta import get_qual, iter_fasta_qual, write_fasta_qual
from jcvi.utils.table import tabulate
from jcvi.apps.softlink import get_abs_path
from jcvi.apps.base import ActionDispatcher, sh, set_grid, debug
from jcvi.assembly.base import CAPATH
debug()


def main():

    actions = (
        ('tracedb', 'convert trace archive files to frg file'),
        ('fasta', 'convert fasta to frg file'),
        ('sff', 'convert 454 reads to frg file'),
        ('fastq', 'convert Illumina reads to frg file'),
        ('script', 'create gatekeeper script file to remove or add mates'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def script(args):
    """
    %prog script bfs_rfs libs

    `bfs_rfs` contains the joined results from Brian's `classifyMates`. We want
    to keep the RFS result (but not in the BFS result) to retain actual MP. Libs
    contain a list of lib iids, use comma to separate, e.g. "9,10,11".
    """
    p = OptionParser(script.__doc__)

    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(p.print_help())

    fsfile, libs = args
    libs = [int(x) for x in libs.split(",")]
    fp = open(fsfile)
    not_found = ("limited", "exhausted")
    counts = defaultdict(int)
    pe, mp = 0, 0
    both, noidea = 0, 0
    total = 0

    for i in libs:
        print "lib iid {0} allfragsunmated 1".format(i)

    for row in fp:
        frgiid, bfs, rfs = row.split()
        bfs = (bfs not in not_found)
        rfs = (rfs not in not_found)
        if bfs and (not rfs):
            pe += 1
        if rfs and (not bfs):
            mp += 1
            frgiid = int(frgiid)
            mateiid = frgiid + 1
            print "frg iid {0} mateiid {1}".format(frgiid, mateiid)
            print "frg iid {0} mateiid {1}".format(mateiid, frgiid)
        if bfs and rfs:
            both += 1
        if (not bfs) and (not rfs):
            noidea += 1
        total += 1

    assert pe + mp + both + noidea == total
    counts[("PE", "N")] = pe
    counts[("MP", "N")] = mp
    counts[("Both", "N")] = both
    counts[("No Idea", "N")] = noidea

    table = tabulate(counts)
    func = lambda a: a * 100. / total
    table = table.withNewColumn("Percentage", callback=func,
            columns=("N",), digits=2)
    print >> sys.stderr, table


def tracedb(args):
    """
    %prog tracedb <xml|lib|frg>

    Run `tracedb-to-frg.pl` within current folder.
    """
    p = OptionParser(tracedb.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    action, = args
    assert action in ("xml", "lib", "frg")

    CMD = "perl {0}tracedb-to-frg.pl".format(CAPATH)
    xmls = glob("xml*")

    if action == "xml":
        for xml in xmls:
            cmd = CMD + " -xml {0}".format(xml)
            sh(cmd, outfile="/dev/null", errfile="/dev/null", background=True)

    elif action == "lib":
        cmd = CMD + " -lib {0}".format(" ".join(xmls))
        sh(cmd)

    elif action == "frg":
        for xml in xmls:
            cmd = CMD + " -frg {0}".format(xml)
            sh(cmd, background=True)


def make_qual(fastafile, defaultqual=20):
    """
    Make a qualfile with default qual value if not available
    """
    assert op.exists(fastafile)

    qualfile = get_qual(fastafile)
    if qualfile is None:
        qualfile = get_qual(fastafile, check=False)
        qualhandle = open(qualfile, "w")

        for rec in iter_fasta_qual(fastafile, None, defaultqual=defaultqual):
            write_fasta_qual(rec, None, qualhandle)

        logging.debug("write qual values to file `{0}`".format(qualfile))
        qualhandle.close()

    return qualfile


def make_matepairs(fastafile):
    """
    Assumes the mates are adjacent sequence records
    """
    assert op.exists(fastafile)

    matefile = fastafile.rsplit(".", 1)[0] + ".mates"
    if op.exists(matefile):
        logging.debug("matepairs file `{0}` found".format(matefile))
    else:
        logging.debug("parsing matepairs from `{0}`".format(fastafile))
        matefw = open(matefile, "w")
        it = SeqIO.parse(fastafile, "fasta")
        for fwd, rev in zip(it, it):
            print >> matefw, "{0}\t{1}".format(fwd.id, rev.id)

        matefw.close()

    return matefile


def add_size_option(p):
    p.add_option("-s", dest="size", default=0, type="int",
            help="insert has mean size of [default: %default] " + \
                 "stddev is assumed to be 25% around mean size")


get_mean_sv = lambda size: (size, size / 4)


def fasta(args):
    """
    %prog fasta fastafile

    Convert reads formatted as FASTA file, and convert to CA frg file. If .qual
    file is found, then use it, otherwise just make a fake qual file. Mates are
    assumed as adjacent sequence records (i.e. /1, /2, /1, /2 ...) unless a
    matefile is given.
    """
    p = OptionParser(fasta.__doc__)
    p.add_option("-m", dest="matefile", default=None,
            help="matepairs file")
    set_grid(p)
    add_size_option(p)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    grid = opts.grid

    fastafile, = args
    libname = op.basename(fastafile).split(".")[0]
    frgfile = fastafile.rsplit(".", 1)[0] + ".frg"

    mated = (opts.size != 0)
    mean, sv = get_mean_sv(opts.size)

    qualfile = make_qual(fastafile)
    if mated:
        if opts.matefile:
            matefile = opts.matefile
            assert op.exists(matefile)
        else:
            matefile = make_matepairs(fastafile)

    cmd = CAPATH + "convert-fasta-to-v2.pl -l {0} -s {1} -q {2} ".\
            format(libname, fastafile, qualfile)
    if mated:
        cmd += "-mean {0} -stddev {1} -m {2} ".format(mean, sv, matefile)

    sh(cmd, grid=grid, outfile=frgfile)


def sff(args):
    """
    %prog sff sffiles

    Convert reads formatted as 454 SFF file, and convert to CA frg file.
    """
    p = OptionParser(sff.__doc__)
    p.add_option("--prefix", dest="prefix", default=None,
            help="output frg filename prefix")
    set_grid(p)
    add_size_option(p)

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(p.print_help())

    grid = opts.grid

    sffiles = args
    plates = [x.split(".")[0].split("_")[-1] for x in sffiles]

    mated = (opts.size != 0)
    mean, sv = get_mean_sv(opts.size)

    if len(plates) > 1:
        plate = plates[0][:-1] + 'X'
    else:
        plate = "_".join(plates)

    if mated:
        libname = "Titan{0}Kb-".format(opts.size / 1000) + plate
    else:
        libname = "TitanFrags-" + plate

    if opts.prefix:
        libname = opts.prefix + "-" + plate

    cmd = CAPATH + "sffToCA -libraryname {0} -output {0} ".format(libname)
    cmd += "-clear 454 -trim chop "
    if mated:
        cmd += "-linker titanium -insertsize {0} {1} ".format(mean, sv)

    cmd += " ".join(sffiles)

    sh(cmd, grid=grid)


def fastq(args):
    """
    %prog fastq fastqfile

    Convert reads formatted as FASTQ file, and convert to CA frg file.
    """
    p = OptionParser(fastq.__doc__)
    p.add_option("--sanger", dest="sanger", default=False, action="store_true",
            help="Are the qv sanger encodings? [default: %default]")
    p.add_option("--outtie", dest="outtie", default=False, action="store_true",
            help="Are these outie reads? [default: %default]")
    p.add_option("--deduplicate", dest="deduplicate", default=False,
            action="store_true", help="set doRemoveDuplicateReads=1 "
            "[default: %default]")
    add_size_option(p)

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(p.print_help())

    fastqfiles = [get_abs_path(x) for x in args]

    mated = (opts.size != 0)
    outtie = opts.outtie
    libname = op.basename(fastqfiles[0]).split(".")[0]
    libname = libname.replace("_1_sequence", "")

    if outtie:
        libname = "IlluminaMP_" + libname
    else:
        libname = "IlluminaPE_" + libname

    if mated:
        libname += "_Mated"
    else:
        if outtie:
            libname = "IlluminaMP_UnMated"
        else:
            libname = "IlluminaPE_UnMated"
    frgfile = libname + ".frg"

    mean, sv = get_mean_sv(opts.size)

    cmd = CAPATH + "fastqToCA -libraryname {0} ".format(libname)
    fastqs = " ".join("-fastq {0}".format(x) for x in fastqfiles)
    if mated:
        assert len(args) == 2, "you need two fastq file for mated library"
        fastqs = "-fastq {0},{1}".format(*fastqfiles)
        cmd += "-insertsize {0} {1} ".format(mean, sv)
    cmd += fastqs

    if opts.sanger:
        cmd += " -type sanger "
    if outtie:
        cmd += " -outtie "

    sh(cmd, outfile=frgfile)
    if opts.deduplicate:
        cmd = "sed -i 's/^doRemoveDuplicateReads.*/doRemoveDuplicateReads=1/' "
        cmd += frgfile
        sh(cmd)


if __name__ == '__main__':
    main()
