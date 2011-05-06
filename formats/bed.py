
"""
Classes to handle the .bed files
"""

import sys
import logging

from jcvi.formats.base import LineFile


class BedLine(object):
    # the Bed format supports more columns. we only need
    # the first 4, but keep the information in 'extra'.
    __slots__ = ("seqid", "start", "end", "accn",
                 "extra", "score", "strand")

    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.seqid = args[0]
        self.start = int(args[1]) + 1
        self.end = int(args[2])
        assert self.start <= self.end, \
                "start={0} end={1}".format(self.start, self.end)
        self.accn = args[3]
        self.extra = self.score = self.strand = None

        if len(args) > 4:
            self.extra = args[4:]
            self.score = self.extra[0]
        if len(args) > 5:
            self.strand = self.extra[1]

    def __str__(self):
        s = "\t".join(str(x) for x in (self.seqid, self.start - 1,
            self.end, self.accn))

        if self.extra:
            s += "\t" + "\t".join(self.extra)
        return s

    def __getitem__(self, key):
        return getattr(self, key)

    def reverse_complement(self, sizes):
        # this function is used in assembly.bundle
        seqid = self.seqid.rstrip('-')
        size = sizes.get_size(seqid)

        if self.seqid[-1] == '-':
            self.seqid = self.seqid[:-1]
        else:
            self.seqid += '-'

        start = size - self.end + 1
        end = size - self.start + 1
        self.start, self.end = start, end
        assert self.start <= self.end, \
                "start={0} end={1}".format(self.start, self.end)

        if self.strand:
            strand = {'+': '-', '-': '+'}[self.strand]
            self.strand = self.extra[1] = strand


class Bed(LineFile):

    def __init__(self, filename=None, key=None):
        super(Bed, self).__init__(filename)

        # the sorting key provides some flexibility in ordering the features
        # for example, user might not like the lexico-order of seqid
        self.key = key or (lambda x: (x.seqid, x.start, x.accn))

        if not filename:
            logging.debug("Initiate bed structure without filename")
            return

        for line in open(filename):
            if line[0] == "#":
                continue
            self.append(BedLine(line))

        self.sort(key=self.key)

    def print_to_file(self, fw=sys.stdout):
        for bedline in self:
            print >>fw, bedline

    @property
    def seqids(self):
        return sorted(set(b.seqid for b in self))

    @property
    def order(self):
        # get the gene order given a Bed object
        return dict((f.accn, (i, f)) for (i, f) in enumerate(self))

    @property
    def simple_bed(self):
        return [(b.seqid, i) for (i, b) in enumerate(self)]

    def sub_bed(self, seqid):
        # get all the beds on one chromosome
        for b in self:
            if b.seqid == seqid:
                yield b
