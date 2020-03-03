#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Classes to handle the .psl files
"""
from __future__ import print_function

import sys
import math
import re

from jcvi.formats.base import LineFile, must_open
from jcvi.apps.base import OptionParser, ActionDispatcher


class PslLine(object):
    def __init__(self, sline):
        args = sline.strip().split()
        self.nargs = len(args)
        self.matches = int(args[0])
        self.misMatches = int(args[1])
        self.repMatches = int(args[2])
        self.nCount = int(args[3])
        self.qNumInsert = int(args[4])
        self.qBaseInsert = int(args[5])
        self.tNumInsert = int(args[6])
        self.tBaseInsert = int(args[7])
        self.qstrand, self.strand = args[8], None
        m = re.match(r"(?P<qs>[\+\-]?)(?P<gs>[\+\-])", self.qstrand)
        if m:
            self.qstrand, self.strand = m.group("qs"), m.group("gs")
        self.qName = args[9]
        self.qSize = int(args[10])
        self.qStart = int(args[11])
        self.qEnd = int(args[12])
        ##        if self.qstrand == "-":
        ##            self.qStart, self.qEnd = self.qSize - self.qEnd, \
        ##                    self.qSize - self.qStart
        self.tName = args[13]
        self.tSize = int(args[14])
        self.tStart = int(args[15])
        self.tEnd = int(args[16])
        self.blockCount = int(args[17])
        self.blockSizes = [int(x) for x in args[18].strip().split(",")[:-1]]
        self.qStarts = [int(x) for x in args[19].strip().split(",")[:-1]]
        self.tStarts = [int(x) for x in args[20].strip().split(",")[:-1]]

    ##        self.tStarts = [self.tSize - int(x) if self.strand == "-" \
    ##                else int(x) for x in args[20].strip().split(',')[:-1]]

    def __str__(self):
        args = [
            self.matches,
            self.misMatches,
            self.repMatches,
            self.nCount,
            self.qNumInsert,
            self.qBaseInsert,
            self.tNumInsert,
            self.tBaseInsert,
            self.strand,
            self.qName,
            self.qSize,
            self.qStart,
            self.qEnd,
            self.tName,
            self.tSize,
            self.tStart,
            self.tEnd,
            self.blockCount,
            self.blockSizes,
            self.qStarts,
            self.tStarts,
        ]

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
        sizeMult = self._sizeMult

        return (
            sizeMult * (self.matches + (self.repMatches >> 1))
            - sizeMult * self.misMatches
            - self.qNumInsert
            - self.tNumInsert
        )

    @property
    def coverage(self):
        return (
            100
            * (self.matches + self.misMatches + self.repMatches + self.nCount)
            / self.qSize
        )

    @property
    def swap(self):
        self.qName, self.qSize, self.tName, self.tSize = (
            self.tName,
            self.tSize,
            self.qName,
            self.qSize,
        )

        self.qStart, self.qEnd, self.tStart, self.tEnd = (
            self.tStart,
            self.tEnd,
            self.qStart,
            self.qEnd,
        )

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
        return (
            (self.tEnd == self.tStarts[last] + 3 * self.blockSizes[last])
            and self.strand == "+"
        ) or (
            (
                self.tStart
                == self.tSize - (self.tStarts[last] + 3 * self.blockSizes[last])
                and self.strand == "-"
            )
        )

    def _milliBad(self, ismRNA=False):
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

        return (
            (
                1000
                * (
                    self.misMatches * sizeMult
                    + insertFactor
                    + round(3 * math.log(1 + sizeDiff))
                )
            )
            / total
            if total != 0
            else 0
        )

    def pct_id(self, simple=None):
        return (
            100.00 - self._milliBad(ismRNA=True) * 0.1
            if not simple
            else 100.00 * self.matches / (self.matches + self.misMatches)
        )
        # else 100.00 * self.score / self.qSize

    def gffline(
        self,
        source="GMAP",
        type="match_part",
        primary_tag="Parent",
        alt_score=None,
        suffix=".match",
        count=0,
    ):

        score = "." if type == "match_part" else "{0:.2f}".format(self.score)

        target = " ".join(str(x) for x in [self.qName, self.qStart, self.qEnd])

        attributes = [
            primary_tag + "=" + self.qName + suffix + str(count),
            "Target=" + target,
        ]
        if primary_tag == "ID":
            attributes.extend(
                [
                    "identity={0:.2f}".format(self.pct_id(simple=alt_score)),
                    "coverage={0:.2f}".format(self.coverage),
                ]
            )
        attrs = ";".join(str(x) for x in attributes)

        line = "\t".join(
            str(x)
            for x in [
                self.tName,
                source,
                type,
                self.tStart,
                self.tEnd,
                score,
                self.strand,
                ".",
                attrs,
            ]
        )
        return line

    @property
    def bed12line(self):
        color = "255,0,0"
        self.blockStarts = ",".join([str(x - self.tStart) for x in self.tStarts])
        line = "\t".join(
            str(x)
            for x in (
                self.tName,
                self.tStart,
                self.tEnd,
                self.qName,
                "{0:.2f}".format(self.pct_id()),
                self.strand,
                self.tStart,
                self.tEnd,
                color,
                self.blockCount,
                ",".join(str(bs) for bs in self.blockSizes),
                self.blockStarts,
            )
        )
        return line


class Psl(LineFile):
    def __init__(self, filename=None):
        super(Psl, self).__init__(filename)

        import re

        self.mCounts = {}  # dict to hold match counts
        if not filename:
            return

        for line in must_open(filename):
            if not re.match(r"\d+", line[0]):
                continue
            self.append(PslLine(line))

    def trackMatches(self, id):
        self.mCounts[id] = self.mCounts.get(id, 0) + 1

    def getMatchCount(self, id):
        return self.mCounts[id]


def main():

    actions = (
        ("gff", "convert psl to gff3 format"),
        ("bed", "convert psl to bed12 format"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed pslfile

    Convert to bed format.
    """
    p = OptionParser(bed.__doc__)
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (pslfile,) = args
    fw = must_open(opts.outfile, "w")

    psl = Psl(pslfile)
    for p in psl:
        print(p.bed12line, file=fw)


def gff(args):
    """
    %prog gff pslfile

    Convert to gff format.
    """
    p = OptionParser(gff.__doc__)
    p.add_option("--source", default="GMAP", help="specify GFF source")
    p.add_option(
        "--type", default="EST_match", help="specify GFF feature type",
    )
    p.add_option("--suffix", default=".match", help="match ID suffix")
    p.add_option(
        "--swap",
        default=False,
        action="store_true",
        help="swap query and target features",
    )
    p.add_option(
        "--simple_score",
        default=False,
        action="store_true",
        help="calculate a simple percent score",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (pslfile,) = args
    fw = must_open(opts.outfile, "w")

    print("##gff-version 3", file=fw)
    psl = Psl(pslfile)
    for p in psl:
        if opts.swap:
            p.swap

        psl.trackMatches(p.qName)
        # switch from 0-origin to 1-origin
        p.qStart += 1
        p.tStart += 1

        print(
            p.gffline(
                source=opts.source,
                type=opts.type,
                suffix=opts.suffix,
                primary_tag="ID",
                alt_score=opts.simple_score,
                count=psl.getMatchCount(p.qName),
            ),
            file=fw,
        )

        # create an empty PslLine() object and load only
        # the targetName, queryName and strand info
        part = PslLine("\t".join(str(x) for x in [0] * p.nargs))
        part.tName, part.qName, part.strand = p.tName, p.qName, p.strand

        nparts = len(p.qStarts)
        for n in range(nparts):
            part.qStart, part.tStart, aLen = (
                p.qStarts[n] + 1,
                p.tStarts[n] + 1,
                p.blockSizes[n],
            )
            part.qEnd = part.qStart + aLen - 1
            part.tEnd = part.tStart + aLen - 1

            if part.strand == "-":
                part.qStart = p.qSize - (p.qStarts[n] + p.blockSizes[n]) + 1
                part.qEnd = p.qSize - p.qStarts[n]

            print(
                part.gffline(
                    source=opts.source,
                    suffix=opts.suffix,
                    count=psl.getMatchCount(part.qName),
                ),
                file=fw,
            )


if __name__ == "__main__":
    main()
