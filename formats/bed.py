
"""
Classes to handle the .bed files
"""

import sys
import shutil
import logging

from optparse import OptionParser

from jcvi.formats.base import LineFile, must_open
from jcvi.utils.cbook import thousands
from jcvi.utils.range import range_union
from jcvi.apps.base import ActionDispatcher, debug, sh
debug()


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
            print >> fw, bedline

    @property
    def sum(self):
        ranges = [(x.seqid, x.start, x.end) for x in self]
        return range_union(ranges)

    @property
    def seqids(self):
        return sorted(set(b.seqid for b in self))

    @property
    def accns(self):
        return sorted(set(b.accn for b in self))

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

    def sub_beds(self):
        # get all the beds on all chromosomes, emitting one at a time
        for bs in self.seqids:
            yield bs, list(self.sub_bed(bs))


def main():

    actions = (
        ('sort', 'sort bed file'),
        ('sum', 'sum the total lengths of the intervals'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def sum(args):
    """
    %prog sum bedfile

    Sum the total lengths of the intervals.
    """
    p = OptionParser(sum.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    bedfile, = args
    bed = Bed(bedfile)
    print >> sys.stderr, "Total seqids: {0}".format(len(bed.seqids))
    print >> sys.stderr, "Total ranges: {0}".format(len(bed))
    print >> sys.stderr, "Total bases: {0} bp".format(thousands(bed.sum))


def sort(args):
    """
    %prog sort bedfile

    Sort bed file to have ascending order of seqid, then start.
    """
    p = OptionParser(sort.__doc__)
    p.add_option("-i", "--inplace", dest="inplace",
            default=False, action="store_true",
            help="Sort bed file in place [default: %default]")
    p.add_option("--fast", dest="fast", default=False, action="store_true",
            help="Use shell `sort` [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    inplace = opts.inplace
    if opts.fast:
        cmd = "sort -k1,1 -k2,2n -k4,4 {0}".format(bedfile)
        if inplace:
            cmd += " -o {0}".format(bedfile)
        sh(cmd)

    else:
        bed = Bed(bedfile)
        sortedbed = bedfile + ".sorted"
        fw = must_open(sortedbed, "w") if inplace else sys.stdout
        bed.print_to_file(fw=fw)

        if inplace:
            shutil.move(sortedbed, bedfile)


if __name__ == '__main__':
    main()
