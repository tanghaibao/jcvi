#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Classes to handle the .psl files
"""

import os
import os.path as op
import sys
import math
import logging

from optparse import OptionParser

from jcvi.formats.base import LineFile, must_open
from jcvi.apps.base import ActionDispatcher, debug, sh, \
        need_update, popen, set_outfile
debug()


class PslLine(object):

    def __init__(self, sline):
        args = sline.strip().split()
        self.nargs = nargs = len(args)

        self.matches = int(args[0])
        self.misMatches = int(args[1])
        self.repMatches = int(args[2])
        self.nCount = int(args[3])
        self.qNumInsert = int(args[4])
        self.qBaseInsert = int(args[5])
        self.tNumInsert = int(args[6])
        self.tBaseInsert = int(args[7])
        self.strand = args[8]
        self.qName = args[9]
        self.qSize = int(args[10])
        self.qStart = int(args[11])
        self.qEnd = int(args[12])
        self.tName = args[13]
        self.tSize = int(args[14])
        self.tStart = int(args[15])
        self.tEnd = int(args[16])
        self.blockCount = int(args[17])
        self.blockSizes = [int(x) for x in args[18].split(',')[0:-1]]
        self.qStarts = [int(x) for x in args[19].split(',')[0:-1]]
        self.tStarts = [int(x) for x in args[20].strip().split(',')[0:-1]]

    def __str__(self):
        args = [self.matches, self.misMatches, self.repMatches, \
                self.nCount, self.qNumInsert, self.qBaseInsert, \
                self.tNumInsert, self.tBaseInsert, self.strand, \
                self.qName, self.qSize, self.qStart, self.qEnd, \
                self.tName, self.tSize, self.tStart, self.tEnd, \
                self.blockCount, self.blockSizes, self.qStarts, \
                self.tStarts]

        s = "\t".join(str(x) for x in args)
        return s

    def __getitem__(self, key):
        return getattr(self, key)

    @property
    def qspan(self):
        return self.qEnd - self.qStart

    @property
    def tspan(self):
        return self.tEnd - self.tStart

    @property
    def score(self):
        return self.matches + (self.repMatches / 2) - self.misMatches - \
                self.qNumInsert - self.tNumInsert

    @property
    def swap(self):
        self.qName, self.qSize, self.tName, self.tSize = \
                self.tName, self.tSize, self.qName, self.qSize

        self.qStart, self.qEnd, self.tStart, self.tEnd = \
                self.tStart, self.tEnd, self.qStart, self.qEnd

        self.qStarts, self.tStarts = self.tStarts, self.qStarts

    @property
    def _sizeMult(self):
        """
        decide the size multiplier based on sequence space (protein/nucleotide)
        """
        return 3 if self._isProtein else 1

    @property
    def _isProtein(self):
        """
        check if blockSizes and scores are in the protein space or not
        """
        last = self.blockCount - 1
        return ((self.tEnd == self.tStarts[last] + 3 * self.blockSizes[last]) \
                and self.strand == "+") or \
                ((self.tSize - (self.tStarts[last] + 3 * self.blockSizes[last])\
                and self.strand == "-"))

    def _milliBad(self, ismRNA=True):
        """
        calculate badness in parts per thousand
        i.e. number of non-identical matches
        """
        sizeMult = self._sizeMult

        qAlnSize, tAlnSize = self.qspan * sizeMult, self.tspan
        alnSize = min(qAlnSize, tAlnSize)
        if alnSize <= 0:
            return 0

        sizeDiff = qAlnSize - tAlnSize
        if sizeDiff < 0:
            sizeDiff = 0 if ismRNA else -sizeDiff

        insertFactor = self.qNumInsert
        if not ismRNA:
            insertFactor += self.tNumInsert

        total = (self.matches + self.repMatches + self.misMatches) * sizeMult

        return (1000 * (round(3 * math.log(1 + sizeDiff) + self.misMatches * sizeMult \
                + insertFactor))) / total if total != 0 else 0

    def pct_id(self, simple=None):
        return 100.00 - self._milliBad(ismRNA=True) * 0.1 if not simple \
                else 100.00 * self.matches / self.qSize

    def gffline(self, source="GMAP", type="match_part", primary_tag="Parent", \
           alt_score=None, suffix=".match", count=0):

        score = "." if type == "match_part" else "{0:.2f}".format(self.pct_id(simple=alt_score))

        target = " ".join(str(x) for x in [self.qName, self.qStart, self.qEnd])
        attributes = ";".join(str(x) for x in [primary_tag + "=" + self.qName + suffix + \
                str(count), "Target=" + target])

        line = "\t".join(str(x) for x in [self.tName, source, type, self.tStart, \
                self.tEnd, score, self.strand, ".", attributes])
        return line


class Psl(LineFile):

    def __init__(self, filename=None):
        super(Psl, self).__init__(filename)

        self.mCounts = {}   # dict to hold match counts
        if not filename:
            return

        for line in must_open(filename):
            if line[0] == "#":
                continue
            self.append(PslLine(line))

    def trackMatches(self, id):
        self.mCounts[id] = self.mCounts.get(id, 0) + 1

    def getMatchCount(self, id):
        return self.mCounts[id]

def main():

    actions = (
        ('gff', 'convert psl to gff3 format'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def gff(args):
    """
    %prog gff pslfile

    Convert to gff format.
    """
    p = OptionParser(gff.__doc__)
    p.add_option("--source", default="GMAP",
                 help="specify GFF source [default: %default]")
    p.add_option("--type", default="EST_match",
                help="specify GFF feature type [default: %default]")
    p.add_option("--suffix", default=".match",
                 help="match ID suffix [default: \"%default\"]")
    p.add_option("--swap", default=False, action="store_true",
                 help="swap query and target features [default: %default]")
    p.add_option("--simple_score", default=False, action="store_true",
                 help="calculate a simple percent score [default: %default]")
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pslfile, = args
    fw = must_open(opts.outfile, "w")

    print >> fw, "##gff-version 3"
    psl = Psl(pslfile)
    for p in psl:
        if opts.swap:
            p.swap

        psl.trackMatches(p.qName)
        # switch from 0-origin to 1-origin
        p.qStart += 1
        p.tStart += 1
        if p.strand == "-":
            p.qStart, p.qEnd = p.qEnd, p.qStart

        print >> fw, p.gffline(source=opts.source, type=opts.type, suffix=opts.suffix, \
                primary_tag="ID", alt_score=opts.simple_score, \
                count=psl.getMatchCount(p.qName))

        # create an empty PslLine() object and load only
        # the targetName, queryName and strand info
        part = PslLine("\t".join(str(x) for x in [0] * p.nargs))
        part.tName, part.qName, part.strand = p.tName, p.qName, p.strand

        nparts = len(p.qStarts)
        for n in xrange(nparts):
            part.qStart, part.tStart, aLen = p.qStarts[n], p.tStarts[n], p.blockSizes[n]
            part.qEnd = part.qStart + aLen
            part.tEnd = part.tStart + aLen

            part.qStart += 1
            part.tStart += 1

            if part.strand == "-":
                part.aLen = p.blockSizes[nparts - 1 - n]
                part.qEnd = p.qStarts[nparts - 1 - n]
                part.qStart = part.qEnd + part.aLen
                part.qEnd += 1

            print >> fw, part.gffline(source=opts.source, suffix=opts.suffix, \
                    count=psl.getMatchCount(part.qName))


if __name__ == '__main__':
    main()
