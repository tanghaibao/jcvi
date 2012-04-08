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
from jcvi.apps.base import ActionDispatcher, debug
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
        ('prepare', 'prepare SOAP config files and run script'),
        ('fillstats', 'build stats on .fill file from GapCloser'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


SOAPRUN="""#!/bin/bash

P=32
S=soap.config
C=SOAPdenovo-63mer
K=47
A=asm${K}

$C pregraph -s $S -K $K -o $A -a 300 -p $P -R
$C contig -g $A -M 3 -R
$C map -p $P -s $S -g $A
$C scaff -F -g $A
GapCloser -t $P -o ${A}.closed.scafSeq -a ${A}.scafSeq -p 31 -b $Si -l 155"""


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
        block += "reverse_seq={0}\n".format(reverse_seq)
        block += "asm_flags={0}\n".format(lib.asm_flags)
        block += "rank={0}\n".format(rank)
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
