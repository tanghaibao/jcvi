#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
TIGR contig format, see spec:

<http://www.cbcb.umd.edu/research/contig_representation.shtml#contig>
"""

import sys
import logging

from jcvi.formats.base import BaseFile, read_block
from jcvi.apps.base import OptionParser, ActionDispatcher


class ReadLine(object):
    def __init__(self, row, contig):
        # '#16(0) [RC] 3046 bases, 00000000 checksum. {3046 1} <1 3046>'
        assert row[0] == "#"
        self.id = row.strip("#").split("(")[0]
        coords = row.split("<")[1].split(">")[0]
        start, end = coords.split()
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        if self.start > self.end:
            self.start, self.end = self.end, self.start

        self.orientation = "-" if "[RC]" in row else "+"

    def __str__(self):
        return self.id

    @property
    def bedline(self):
        return "\t".join(
            str(x)
            for x in (
                self.contig,
                self.start - 1,
                self.end,
                self.id,
                "0",
                self.orientation,
            )
        )

    __repr__ = __str__


class ContigLine(object):
    def __init__(self, row):
        # '##1 6 8914 bases, 00000000 checksum.'
        assert row[:2] == "##"
        self.id = row.strip("#").split()[0]
        self.reads = []

    def __str__(self):
        return ":".join((self.id, str(self.reads)))

    __repr__ = __str__


class ContigFile(BaseFile):
    def __init__(self, filename):
        super(ContigFile, self).__init__(filename)
        self.fp = open(filename)

    def iter_records(self):
        c = None
        for a, b in read_block(self.fp, "#"):
            if a[:2] == "##":
                if c:
                    yield c
                c = ContigLine(a)
            else:
                c.reads.append(ReadLine(a, c.id))
        if c:  # last one
            yield c


def main():

    actions = (
        ("bed", "convert read membership to bed format"),
        ("frombed", "convert read placement to contig format"),
    )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def frombed(args):
    """
    %prog frombed bedfile contigfasta readfasta

    Convert read placement to contig format. This is useful before running BAMBUS.
    """
    from jcvi.formats.fasta import Fasta
    from jcvi.formats.bed import Bed
    from jcvi.utils.cbook import fill

    p = OptionParser(frombed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, contigfasta, readfasta = args
    prefix = bedfile.rsplit(".", 1)[0]
    contigfile = prefix + ".contig"
    idsfile = prefix + ".ids"

    contigfasta = Fasta(contigfasta)
    readfasta = Fasta(readfasta)

    bed = Bed(bedfile)
    checksum = "00000000 checksum."
    fw_ids = open(idsfile, "w")
    fw = open(contigfile, "w")

    for ctg, reads in bed.sub_beds():
        ctgseq = contigfasta[ctg]
        ctgline = "##{0} {1} {2} bases, {3}".format(
            ctg, len(reads), len(ctgseq), checksum
        )

        print(ctg, file=fw_ids)
        print(ctgline, file=fw)
        print(fill(ctgseq.seq), file=fw)

        for b in reads:
            read = b.accn
            strand = b.strand
            readseq = readfasta[read]
            rc = " [RC]" if strand == "-" else ""
            readlen = len(readseq)
            rstart, rend = 1, readlen
            if strand == "-":
                rstart, rend = rend, rstart

            readrange = "{{{0} {1}}}".format(rstart, rend)
            conrange = "<{0} {1}>".format(b.start, b.end)
            readline = "#{0}(0){1} {2} bases, {3} {4} {5}".format(
                read, rc, readlen, checksum, readrange, conrange
            )
            print(readline, file=fw)
            print(fill(readseq.seq), file=fw)

    logging.debug("Mapped contigs written to `{0}`.".format(contigfile))
    logging.debug("Contig IDs written to `{0}`.".format(idsfile))


def bed(args):
    """
    %prog bed contigfile

    Prints out the contigs and their associated reads.
    """
    p = OptionParser(main.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (contigfile,) = args
    bedfile = contigfile.rsplit(".", 1)[0] + ".bed"
    fw = open(bedfile, "w")
    c = ContigFile(contigfile)

    for rec in c.iter_records():
        for r in rec.reads:
            print(r.bedline, file=fw)

    logging.debug("File written to `%s`.", bedfile)

    return bedfile


if __name__ == "__main__":
    main()
