#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Script to write and assist SOAPdenovo assembly.
"""

import os.path as op
import sys
import logging

from itertools import groupby
from optparse import OptionParser

from jcvi.assembly.base import FastqNamings, Library
from jcvi.apps.base import ActionDispatcher, debug, need_update, sh
debug()


class FillLine (object):

    def __init__(self, row):
        args = row.split()
        self.start = int(args[0])
        self.end = int(args[1])
        self.leftextend = int(args[2])
        self.rightextend = int(args[3])
        self.closed = (int(args[4]) == 1)
        self.extendlength = int(args[5])
        self.before = int(args[6])
        self.after = int(args[7])
        # Convert from unsigned to signed
        # <http://stackoverflow.com/questions/1375897/how-to-get-the-signed-integer-value-of-a-long-in-python>
        if (self.after & 0x80000000):
            self.after += -0x100000000

    @property
    def delta(self):
        return self.after - self.before


def main():

    actions = (
        ('clean', 'clean and dedup paired FASTQ files'),
        ('correct', 'correct reads using ErrorCorrection'),
        ('prepare', 'prepare SOAP config files and run script'),
        ('fillstats', 'build stats on .fill file from GapCloser'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


SOAPRUN="""#!/bin/bash

P=32
S=soap.config
C=SOAPdenovo-63mer
K=29
A=asm${K}

$C pregraph -s $S -d 1 -p $P -K $K -o $A -p $P
$C contig -g $A -M 3
$C map -s $S -g $A -p $P
$C scaff -g $A -b 1.2 -F -p $P
GapCloser -a ${A}.scafSeq -b $S -l 155 -o ${A}.closed.scafSeq -p 29 -t $P"""


def correct(args):
    """
    %prog correct *.dedup

    Correct reads using ErrorCorrection. Only PE will be used to build the K-mer
    table, but both PE and MP will be corrected. Final command needs to be run:

    $ parallel ErrorCorrection correct output.freq.gz output.sfreq.gz {} ::: xa?
    """
    p = OptionParser(correct.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 3:
        sys.exit(not p.print_help())

    lstfile = "reads2cor.lst"
    fw = open(lstfile, "w")
    print >> fw, "\n".join(x for x in args if x[:2] == "PE")
    fw.close()

    if need_update(args, "output.freq.gz"):
        cmd = "ErrorCorrection kmerfreq reads2cor.lst"
        sh(cmd)

    fw = open(lstfile, "w")
    print >> fw, "\n".join(args)
    fw.close()

    sh("split -l2 {0}".format(lstfile))


def clean(args):
    """
    %prog clean 1.fastq 2.fastq insertsize

    Clean and dedup paired FASTQ files.
    """
    from jcvi.formats.fastq import guessoffset, convert

    p = OptionParser(clean.__doc__)
    p.add_option("-a", default=0, type="int",
                 help="Trim length at 5' end [default: %default]")
    p.add_option("-b", default=0, type="int",
                 help="Trim length at 5' end [default: %default]")
    p.add_option("--nofragsdedup", default=False, action="store_true",
                 help="Don't deduplicate the fragment reads [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    p1, p2, size = args
    size = int(size)
    pf = p1.split(".")[0]

    offset = guessoffset([p1])
    p1_q64 = p1.replace(".gz", "") + ".q64.gz"
    p2_q64 = p2.replace(".gz", "") + ".q64.gz"
    a, b = opts.a, opts.b

    if offset == 33 and need_update([p1, p2], [p1_q64, p2_q64]):
        logging.debug("Converting offset from 33 to 64 ...")
        convert([p1, p1_q64, "-Q", "sanger", "-q", "illumina"])
        convert([p2, p2_q64, "-Q", "sanger", "-q", "illumina"])
        p1, p2 = p1_q64, p2_q64

    p1_clean = p1 + ".clean"
    p2_clean = p2 + ".clean"
    if need_update([p1, p2], [p1_clean, p2_clean]):
        logging.debug("Running low quality filtering ...")
        cmd  = "filter_data_gz -y -z -w 10 -B 40"
        cmd += " -l {0} -a {1} -b {2} -c {1} -d {2}".format(size, a, b, a, b)
        cmd += " {0} {1} {2}.clean.stat {3} {4}".\
                    format(p1, p2, pf, p1_clean, p2_clean)
        sh(cmd)

    if opts.nofragsdedup:
        return

    p1, p2 = p1_clean, p2_clean
    p1_dedup = p1 + ".dedup"
    p2_dedup = p2 + ".dedup"
    if need_update([p1, p2], [p1_dedup, p2_dedup]):
        logging.debug("Running duplicate filtering ...")
        cmd  = "duplication {0} {1} {2} {3} {4}.dedup.stat".\
                    format(p1, p2, p1_dedup, p2_dedup, pf)
        sh(cmd)


def fillstats(args):
    """
    %prog fillstats genome.fill

    Build stats on .fill file from GapCloser.
    """
    from jcvi.utils.cbook import SummaryStats, percentage, thousands

    p = OptionParser(fillstats.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fillfile, = args
    fp = open(fillfile)
    scaffolds = 0
    gaps = []
    for row in fp:
        if row[0] == ">":
            scaffolds += 1
            continue
        fl = FillLine(row)
        gaps.append(fl)

    print >> sys.stderr, "{0} scaffolds in total".format(scaffolds)

    closed = [x for x in gaps if x.closed]
    closedbp = sum(x.before for x in closed)
    notClosed = [x for x in gaps if not x.closed]
    notClosedbp = sum(x.before for x in notClosed)

    totalgaps = len(closed) + len(notClosed)
    totalbp = closedbp + notClosedbp

    print >> sys.stderr, "Closed gaps: {0} size: {1} bp".\
                        format(percentage(len(closed), totalgaps), thousands(closedbp))
    ss = SummaryStats([x.after for x in closed])
    print >> sys.stderr, ss

    ss = SummaryStats([x.delta for x in closed])
    print >> sys.stderr, "Delta:", ss

    print >> sys.stderr, "Remaining gaps: {0} size: {1} bp".\
                        format(percentage(len(notClosed), totalgaps), thousands(notClosedbp))
    ss = SummaryStats([x.after for x in notClosed])
    print >> sys.stderr, ss


def prepare(args):
    """
    %prog prepare *.fastq

    Scan input fastq files (see below) and write SOAP config files based
    on inputfiles.
    """
    from jcvi.utils.iter import grouper
    from jcvi.formats.base import check_exists

    p = OptionParser(prepare.__doc__ + FastqNamings)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fnames = args
    for x in fnames:
        assert op.exists(x), "File `{0}` not found.".format(x)

    cfgfile = "soap.config"
    fw = open(cfgfile, "w")

    library_name = lambda x: "-".join(\
                op.basename(x).split(".")[0].split("-")[:2])
    libs = [(Library(x), sorted(fs)) for x, fs in \
                groupby(fnames, key=library_name)]

    libs.sort(key=lambda x: x[0].size)
    rank = 0
    singletons = []
    for lib, fs in libs:
        size = lib.size
        if size == 0:
            singletons = fs
            continue

        rank += 1
        block = "[LIB]\n"
        block += "avg_ins={0}\n".format(size)
        f = fs[0]
        reverse_seq = 0 if ".corr." in f else lib.reverse_seq
        pair_num_cutoff = 5 if lib.reverse_seq else 3
        block += "reverse_seq={0}\n".format(reverse_seq)
        block += "asm_flags={0}\n".format(lib.asm_flags)
        block += "rank={0}\n".format(rank)
        block += "pair_num_cutoff={0}\n".format(pair_num_cutoff)
        if lib.reverse_seq:
            block += "map_len=35\n"

        if singletons:
            fs += singletons
            singletons = []

        for f in fs:
            if ".1." in f:
                tag = "q1"
            elif ".2." in f:
                tag = "q2"
            else:
                tag = "q"
            block += "{0}={1}\n".format(tag, f)
        print >> sys.stderr, block
        print >> fw, block

    runfile = "run.sh"
    if check_exists(runfile):
        fw = open(runfile, "w")
        print >> fw, SOAPRUN
        logging.debug("Run script written to `{0}`.".format(runfile))


if __name__ == '__main__':
    main()
