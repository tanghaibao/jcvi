"""
Classes to handle the .bed files
"""

import os.path as op
import sys
import shutil
import logging

from itertools import groupby
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
        self.nullkey = lambda x: (x.seqid, x.start, x.accn)
        self.key = key or self.nullkey

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

        self.sort(key=self.nullkey)
        # get all the beds on all chromosomes, emitting one at a time
        for bs, sb in groupby(self, key=lambda x: x.seqid):
            yield bs, list(sb)


def main():

    actions = (
        ('sort', 'sort bed file'),
        ('sum', 'sum the total lengths of the intervals'),
        ('pairs', 'estimate insert size between paired reads from bedfile'),
        ('mates', 'print paired reads from bedfile'),
        ('sizes', 'infer the sizes for each seqid'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def sizes(args):
    """
    %prog sizes bedfile

    Infer the sizes for each seqid. Useful before dot plots.
    """
    p = OptionParser(sizes.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    assert op.exists(bedfile)

    sizesfile = bedfile.rsplit(".", 1)[0] + ".sizes"

    fw = must_open(sizesfile, "w", checkexists=True, skipcheck=True)
    if fw:
        b = Bed(bedfile)
        for s, sbeds in b.sub_beds():
            print >> fw, "{0}\t{1}".format(\
                         s, max(x.end for x in sbeds))
        logging.debug("Sizes file written to `{0}`.".format(sizesfile))

    return sizesfile


def pairs(args):
    """
    See __doc__ for set_options_pairs().
    """
    from jcvi.formats.blast import report_pairs, set_options_pairs
    p = set_options_pairs()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args

    basename = bedfile.split(".")[0]
    insertsfile = ".".join((basename, "inserts"))

    sortedbedfile = op.basename(bedfile).rsplit(".", 1)[0] + ".sorted.bed"
    if not op.exists(sortedbedfile):
        bedfile = sort([bedfile, "--accn"])
    else:
        bedfile = sortedbedfile

    fp = open(bedfile)
    data = [BedLine(row) for i, row in enumerate(fp) if i < opts.nrows]

    ascii = not opts.pdf
    return bedfile, report_pairs(data, opts.cutoff, opts.mateorientation,
           pairsfile=opts.pairsfile, insertsfile=insertsfile,
           rclip=opts.rclip, ascii=ascii, bins=opts.bins,
           distmode=opts.distmode)


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

    Sort bed file to have ascending order of seqid, then start. It uses the
    `sort` command.
    """
    p = OptionParser(sort.__doc__)
    p.add_option("-i", "--inplace", dest="inplace",
            default=False, action="store_true",
            help="Sort bed file in place [default: %default]")
    p.add_option("--accn", default=False, action="store_true",
            help="Sort based on the accessions [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    inplace = opts.inplace

    sortedbed = op.basename(bedfile).rsplit(".", 1)[0] + ".sorted.bed"
    if inplace:
        sortedbed = bedfile

    sortopt = "-k1,1 -k2,2n -k4,4" if not opts.accn else \
              "-k4,4 -k1,1 -k2,2n"
    cmd = "sort {0} {1}".format(sortopt, bedfile)
    cmd += " -o {0}".format(sortedbed)
    sh(cmd)

    return sortedbed


def mates(args):
    """
    %prog mates bedfile

    Generate the mates file by inferring from the names.
    """
    p = OptionParser(mates.__doc__)
    p.add_option("--lib", default=None,
            help="Output library information along with pairs [default: %default]")
    p.add_option("--nointra", default=False, action="store_true",
            help="Remove mates that are intra-scaffold [default: %default]")
    p.add_option("--prefix", default=False, action="store_true",
            help="Only keep links between IDs with same prefix [default: %default]")
    p.add_option("--rclip", default=1, type="int",
            help="Pair ID is derived from rstrip N chars [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    rclip = opts.rclip
    lib = opts.lib

    key = (lambda x: x.accn[:-rclip]) if rclip else (lambda x: x.accn)
    bed = Bed(bedfile, key=key)

    pf = bedfile.rsplit(".", 1)[0]
    matesfile = pf + ".mates"
    fw = open(matesfile, "w")
    if lib:
        bedfile, (meandist, stdev, p0, p1, p2) = pairs([bedfile, \
                "--mateorientation=+-"])
        sv = int(1.97 * stdev)
        print >> fw, "\t".join(str(x) for x in \
                ("library", lib, meandist -  sv, meandist + sv))

    num_fragments = num_pairs = 0
    matesbedfile = matesfile + ".bed"
    fwm = open(matesbedfile, "w")
    for pe, lines in groupby(bed, key=key):
        lines = list(lines)
        if len(lines) != 2:
            num_fragments += len(lines)
            continue

        a, b = lines

        if opts.nointra and a.seqid == b.seqid:
            continue

        # Use --prefix to limit the links between seqids with the same prefix
        # For example, contigs of the same BAC, mth2-23j10_001, mth-23j10_002
        if opts.prefix:
            aprefix = a.seqid.split("_")[0]
            bprefix = b.seqid.split("_")[0]
            if aprefix != bprefix:
                continue

        num_pairs += 1
        pair = [a.accn, b.accn]
        if lib:
            pair.append(lib)
        print >> fw, "\t".join(pair)

        print >> fwm, a
        print >> fwm, b

    logging.debug("Discard {0} frags and write {1} pairs to `{2}` and `{3}`.".\
            format(num_fragments, num_pairs, matesfile, matesbedfile))

    fw.close()
    fwm.close()

    return matesfile, matesbedfile


if __name__ == '__main__':
    main()
