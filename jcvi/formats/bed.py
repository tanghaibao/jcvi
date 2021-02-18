"""
Classes to handle the .bed files
"""
import os
import os.path as op
import sys
import math
import logging
import numpy as np

from collections import defaultdict, OrderedDict
from itertools import groupby

from more_itertools import pairwise
from natsort import natsorted, natsort_key

from jcvi.formats.base import DictFile, LineFile, must_open, is_number, get_number
from jcvi.formats.sizes import Sizes
from jcvi.utils.cbook import SummaryStats, thousands, percentage
from jcvi.utils.grouper import Grouper
from jcvi.utils.range import (
    Range,
    range_union,
    range_chain,
    range_distance,
    range_intersect,
)
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, need_update, popen


class BedLine(object):
    # the Bed format supports more columns. we only need
    # the first 4, but keep the information in 'extra'.
    __slots__ = (
        "seqid",
        "start",
        "end",
        "accn",
        "extra",
        "score",
        "strand",
        "args",
        "nargs",
    )

    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.nargs = nargs = len(args)
        self.seqid = args[0]
        self.start = int(args[1]) + 1
        self.end = int(args[2])
        assert self.start <= self.end, "start={0} end={1}".format(self.start, self.end)
        self.extra = self.accn = self.score = self.strand = None

        if nargs > 3:
            self.accn = args[3]
        if nargs > 4:
            self.score = args[4]
        if nargs > 5:
            self.strand = args[5]
        if nargs > 6:
            self.extra = args[6:]

        self.args = args

    def __str__(self):
        args = [self.seqid, self.start - 1, self.end]
        if self.accn is not None:
            args += [self.accn]
        if self.score is not None:
            args += [self.score]
        if self.strand is not None:
            args += [self.strand]
        if self.extra is not None:
            args += self.extra

        s = "\t".join(str(x) for x in args)
        return s

    __repr__ = __str__

    def __getitem__(self, key):
        return getattr(self, key)

    @property
    def span(self):
        return self.end - self.start + 1

    @property
    def range(self):
        strand = self.strand or "+"
        return (self.seqid, self.start, self.end, strand)

    @property
    def tag(self):
        return "{0}:{1}-{2}".format(self.seqid, self.start, self.end)

    def reverse_complement(self, sizes):
        size = sizes.get_size(self.seqid)

        start = size - self.end + 1
        end = size - self.start + 1
        self.start, self.end = start, end
        assert self.start <= self.end, "start={0} end={1}".format(self.start, self.end)

        if self.strand:
            strand = {"+": "-", "-": "+"}[self.strand]

    def gffline(self, type="match", source="default"):
        score = (
            "."
            if not self.score or (self.score and not is_number(self.score))
            else self.score
        )
        strand = "." if not self.strand else self.strand
        row = "\t".join(
            (
                self.seqid,
                source,
                type,
                str(self.start),
                str(self.end),
                score,
                strand,
                ".",
                "ID=" + self.accn,
            )
        )
        return row


class Bed(LineFile):
    def __init__(self, filename=None, key=None, sorted=True, juncs=False, include=None):
        super(Bed, self).__init__(filename)

        # the sorting key provides some flexibility in ordering the features
        # for example, user might not like the lexico-order of seqid
        self.nullkey = lambda x: (natsort_key(x.seqid), x.start, x.accn)
        self.key = key or self.nullkey

        if not filename:
            return

        for line in must_open(filename):
            if line[0] == "#" or (juncs and line.startswith("track name")):
                continue
            b = BedLine(line)
            if include and b.accn not in include:
                continue
            self.append(b)

        if sorted:
            self.sort(key=self.key)

    def add(self, row):
        self.append(BedLine(row))

    def print_to_file(self, filename="stdout", sorted=False):
        if sorted:
            self.sort(key=self.key)

        fw = must_open(filename, "w")
        for b in self:
            if b.start < 1:
                logging.error("Start < 1. Reset start for `{0}`.".format(b.accn))
                b.start = 1
            print(b, file=fw)
        fw.close()

    def sum(self, seqid=None, unique=True):
        return bed_sum(self, seqid=seqid, unique=unique)

    @property
    def seqids(self):
        return natsorted(set(b.seqid for b in self))

    @property
    def accns(self):
        return natsorted(set(b.accn for b in self))

    @property
    def order(self):
        # get the gene order given a Bed object
        return dict((f.accn, (i, f)) for (i, f) in enumerate(self))

    @property
    def order_in_chr(self):
        # get the gene order on a particular seqid
        res = {}
        self.sort(key=self.nullkey)
        for seqid, beds in groupby(self, key=lambda x: x.seqid):
            for i, f in enumerate(beds):
                res[f.accn] = (seqid, i, f)
        return res

    @property
    def bp_in_chr(self):
        # get the bp position on a particular seqid
        res = {}
        self.sort(key=self.nullkey)
        for seqid, beds in groupby(self, key=lambda x: x.seqid):
            for i, f in enumerate(beds):
                res[f.accn] = (seqid, (f.start + f.end) / 2, f)
        return res

    @property
    def max_bp_in_chr(self):
        # Get the maximum bp position on a particular seqid
        res = OrderedDict()
        self.sort(key=self.nullkey)
        for seqid, beds in groupby(self, key=lambda x: x.seqid):
            res[seqid] = max(x.end for x in beds)
        return res

    @property
    def simple_bed(self):
        return [(b.seqid, i) for (i, b) in enumerate(self)]

    @property
    def links(self):
        r = []
        for s, sb in self.sub_beds():
            for a, b in pairwise(sb):
                r.append(((a.accn, a.strand), (b.accn, b.strand)))
        return r

    def extract(self, seqid, start, end):
        # get all features within certain range
        for b in self:
            if b.seqid != seqid:
                continue
            if b.start < start or b.end > end:
                continue
            yield b

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

    def get_breaks(self):
        # get chromosome break positions
        simple_bed = self.simple_bed
        for seqid, ranks in groupby(simple_bed, key=lambda x: x[0]):
            ranks = list(ranks)
            # chromosome, extent of the chromosome
            yield seqid, ranks[0][1], ranks[-1][1]


class BedpeLine(object):
    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.seqid1 = args[0]
        self.start1 = int(args[1]) + 1
        self.end1 = int(args[2])
        self.seqid2 = args[3]
        self.start2 = int(args[4]) + 1
        self.end2 = int(args[5])
        self.accn = args[6]
        self.score = args[7]
        self.strand1 = args[8]
        self.strand2 = args[9]
        self.isdup = False

    @property
    def innerdist(self):
        if self.seqid1 != self.seqid2:
            return -1
        return abs(self.start2 - self.end1)

    @property
    def outerdist(self):
        if self.seqid1 != self.seqid2:
            return -1
        return abs(self.end2 - self.start1)

    @property
    def is_innie(self):
        return (self.strand1, self.strand2) == ("+", "-")

    def rc(self):
        self.strand1 = "+" if self.strand1 == "-" else "-"
        self.strand2 = "+" if self.strand2 == "-" else "-"

    def _extend(self, rlen, size, start, end, strand):
        if strand == "+":
            end = start + rlen - 1
            if end > size:
                end = size
                start = end - rlen + 1
        else:
            start = end - rlen + 1
            if start < 1:
                start = 1
                end = start + rlen - 1
        return start, end, strand

    def extend(self, rlen, size):
        self.start1, self.end1, self.strand1 = self._extend(
            rlen, size, self.start1, self.end1, self.strand1
        )
        self.start2, self.end2, self.strand2 = self._extend(
            rlen, size, self.start2, self.end2, self.strand2
        )

    def __str__(self):
        args = (
            self.seqid1,
            self.start1 - 1,
            self.end1,
            self.seqid2,
            self.start2 - 1,
            self.end2,
            self.accn,
            self.score,
            self.strand1,
            self.strand2,
        )
        return "\t".join(str(x) for x in args)

    @property
    def bedline(self):
        assert self.seqid1 == self.seqid2
        assert self.start1 <= self.end2
        args = (self.seqid1, self.start1 - 1, self.end2, self.accn)
        return "\t".join(str(x) for x in args)


class BedEvaluate(object):
    def __init__(self, TPbed, FPbed, FNbed, TNbed):

        self.TP = Bed(TPbed).sum(unique=True)
        self.FP = Bed(FPbed).sum(unique=True)
        self.FN = Bed(FNbed).sum(unique=True)
        self.TN = Bed(TNbed).sum(unique=True)

    def __str__(self):
        from jcvi.utils.table import tabulate

        table = {}
        table[("Prediction-True", "Reality-True")] = self.TP
        table[("Prediction-True", "Reality-False")] = self.FP
        table[("Prediction-False", "Reality-True")] = self.FN
        table[("Prediction-False", "Reality-False")] = self.TN
        msg = str(tabulate(table))

        msg += "\nSensitivity [TP / (TP + FN)]: {0:.1f} %\n".format(
            self.sensitivity * 100
        )
        msg += "Specificity [TP / (TP + FP)]: {0:.1f} %\n".format(
            self.specificity * 100
        )
        msg += "Accuracy [(TP + TN) / (TP + FP + FN + TN)]: {0:.1f} %".format(
            self.accuracy * 100
        )
        return msg

    @property
    def sensitivity(self):
        if self.TP + self.FN == 0:
            return 0
        return self.TP * 1.0 / (self.TP + self.FN)

    @property
    def specificity(self):
        if self.TP + self.FP == 0:
            return 0
        return self.TP * 1.0 / (self.TP + self.FP)

    @property
    def accuracy(self):
        if self.TP + self.FP + self.FN + self.TN == 0:
            return 0
        return (self.TP + self.TN) * 1.0 / (self.TP + self.FP + self.FN + self.TN)

    @property
    def score(self):
        return "|".join(
            (
                "{0:.3f}".format(x)
                for x in (self.sensitivity, self.specificity, self.accuracy)
            )
        )


class BedSummary(object):
    def __init__(self, bed):
        mspans = [(x.span, x.accn) for x in bed]
        spans, accns = zip(*mspans)
        self.mspans = mspans
        self.stats = SummaryStats(spans)
        self.nseqids = len(set(x.seqid for x in bed))
        self.nfeats = len(bed)
        self.total_bases = bed_sum(bed, unique=False)
        self.unique_bases = bed_sum(bed)
        self.coverage = self.total_bases * 1.0 / self.unique_bases

    def report(self):
        print("Total seqids: {0}".format(self.nseqids), file=sys.stderr)
        print("Total ranges: {0}".format(self.nfeats), file=sys.stderr)
        print(
            "Total unique bases: {0} bp".format(thousands(self.unique_bases)),
            file=sys.stderr,
        )
        print(
            "Total bases: {0} bp".format(thousands(self.total_bases)), file=sys.stderr
        )
        print("Estimated coverage: {0:.1f}x".format(self.coverage), file=sys.stderr)
        print(self.stats, file=sys.stderr)
        maxspan, maxaccn = max(self.mspans)
        minspan, minaccn = min(self.mspans)
        print("Longest: {0} ({1})".format(maxaccn, maxspan), file=sys.stderr)
        print("Shortest: {0} ({1})".format(minaccn, minspan), file=sys.stderr)

    def __str__(self):
        return "\t".join(str(x) for x in (self.nfeats, self.unique_bases))


def bed_sum(beds, seqid=None, unique=True):
    if seqid:
        ranges = [(x.seqid, x.start, x.end) for x in beds if x.seqid == seqid]
    else:
        ranges = [(x.seqid, x.start, x.end) for x in beds]

    unique_sum = range_union(ranges)
    raw_sum = sum(x.span for x in beds)
    return unique_sum if unique else raw_sum


def main():

    actions = (
        ("depth", "calculate average depth per feature using coverageBed"),
        ("mergebydepth", "returns union of features beyond certain depth"),
        ("sort", "sort bed file"),
        ("merge", "merge bed files"),
        ("index", "index bed file using tabix"),
        ("bins", "bin bed lengths into each window"),
        ("summary", "summarize the lengths of the intervals"),
        ("evaluate", "make truth table and calculate sensitivity and specificity"),
        ("pile", "find the ids that intersect"),
        ("pairs", "estimate insert size between paired reads from bedfile"),
        ("mates", "print paired reads from bedfile"),
        ("sizes", "infer the sizes for each seqid"),
        ("uniq", "remove overlapping features with higher scores"),
        ("longest", "select longest feature within overlapping piles"),
        ("bedpe", "convert to bedpe format"),
        ("distance", "calculate distance between bed features"),
        ("sample", "sample bed file and remove high-coverage regions"),
        ("refine", "refine bed file using a second bed file"),
        ("flanking", "get n flanking features for a given position"),
        ("some", "get a subset of bed features given a list"),
        ("fix", "fix non-standard bed files"),
        ("filter", "filter bedfile to retain records between size range"),
        ("filterbedgraph", "filter bedgraph to extract unique regions"),
        ("random", "extract a random subset of features"),
        ("juncs", "trim junctions.bed overhang to get intron, merge multiple beds"),
        ("seqids", "print out all seqids on one line"),
        ("alignextend", "alignextend based on BEDPE and FASTA ref"),
        ("clr", "extract clear range based on BEDPE"),
        ("chain", "chain bed segments together"),
        ("density", "calculates density of features per seqid"),
        ("tiling", "compute the minimum tiling path"),
        ("format", "reformat BED file"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def format(args):
    """
    %prog format input.bed

    Re-format BED file, e.g. switch sequence ids.
    """
    p = OptionParser(format.__doc__)
    p.add_option("--prefix", help="Add prefix to name column (4th)")
    p.add_option(
        "--switch", help="Switch seqids based on two-column file"
    )
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    switch = DictFile(opts.switch, delimiter="\t") if opts.switch else None
    prefix = opts.prefix
    bed = Bed(bedfile)
    with must_open(opts.outfile, "w") as fw:
        for b in bed:
            if prefix:
                b.accn = prefix + b.accn
            if switch and b.seqid in switch:
                b.seqid = switch[b.seqid]
            print(b, file=fw)


def filterbedgraph(args):
    """
    %prog filterbedgraph a.bedgraph 1

    Filter the bedGraph, typically from the gem-mappability pipeline. Unique
    regions are 1, two copies .5, etc.
    """
    p = OptionParser(filterbedgraph.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedgraphfile, cutoff = args
    c = float(cutoff)
    fp = open(bedgraphfile)
    pf = bedgraphfile.rsplit(".", 1)[0]
    filteredbed = pf + ".filtered-{}.bed".format(cutoff)
    fw = open(filteredbed, "w")
    nfiltered = ntotal = 0
    for row in fp:
        b = BedLine(row)
        ntotal += 1
        if float(b.accn) >= c:
            print(b, file=fw)
            nfiltered += 1
    fw.close()
    logging.debug(
        "A total of {} intervals (score >= {}) written to `{}`".format(
            percentage(nfiltered, ntotal), cutoff, filteredbed
        )
    )

    mergeBed(filteredbed, sorted=True, delim=None)


def tiling(args):
    """
    %prog tiling bedfile

    Compute minimum tiling path using as few clones as possible. Implemented
    with dynamic programming. Greedy algorithm may also work according a
    stackoverflow source.
    """
    p = OptionParser(tiling.__doc__)
    p.add_option(
        "--overlap",
        default=3000,
        type="int",
        help="Minimum amount of overlaps required",
    )
    p.set_verbose()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    ov = opts.overlap

    bed = Bed(bedfile)
    inf = len(bed)
    selected = Bed()
    for seqid, sbed in bed.sub_beds():
        g = Grouper()
        current = sbed[0]
        # Partition connected features
        for a in sbed:
            g.join(a)
            # requires a real overlap
            if a.start < current.end - ov:
                g.join(a, current)
            if a.end > current.end:
                current = a

        # Process per partition
        for gbed in g:
            end = max(x.end for x in gbed)
            gbed.sort(key=lambda x: (x.start, -x.end))
            entries = len(gbed)
            counts = [inf] * entries
            counts[0] = 1
            traceback = [-1] * entries
            for i, a in enumerate(gbed):
                for j in range(i + 1, entries):
                    b = gbed[j]
                    if b.start >= a.end - ov:
                        break
                    # Two ranges overlap!
                    if counts[i] + 1 < counts[j]:
                        counts[j] = counts[i] + 1
                        traceback[j] = i
            endi = [i for i, a in enumerate(gbed) if a.end == end]
            last = min((traceback[i], i) for i in endi)[1]
            chain = []
            while last != -1:
                chain.append(last)
                last = traceback[last]
            chain = chain[::-1]
            selected.extend([gbed[x] for x in chain])

            if opts.verbose:
                print(counts)
                print(traceback)
                print(chain)
                print("\n".join(str(x) for x in gbed))
                print("*" * 30)
                print("\n".join(str(gbed[x]) for x in chain))
                print()

    tilingbedfile = bedfile.rsplit(".", 1)[0] + ".tiling.bed"
    selected.print_to_file(filename=tilingbedfile, sorted=True)
    logging.debug(
        "A total of {} tiling features written to `{}`".format(
            len(selected), tilingbedfile
        )
    )


def chain(args):
    """
    %prog chain bedfile

    Chain BED segments together.
    """
    p = OptionParser(chain.__doc__)
    p.add_option("--dist", default=100000, help="Chaining distance")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    cmd = "sort -k4,4 -k1,1 -k2,2n -k3,3n {0} -o {0}".format(bedfile)
    sh(cmd)
    bed = Bed(bedfile, sorted=False)
    newbed = Bed()
    for accn, bb in groupby(bed, key=lambda x: x.accn):
        bb = list(bb)
        g = Grouper()
        for a in bb:
            g.join(a)
        for a, b in pairwise(bb):
            if a.seqid == b.seqid and b.start - a.end < opts.dist:
                g.join(a, b)
        data = []
        for p in g:
            seqid = p[0].seqid
            start = min(x.start for x in p)
            end = max(x.end for x in p)
            score = sum(x.span for x in p)
            data.append((seqid, start - 1, end, accn, score))

        d = max(data, key=lambda x: x[-1])
        newbed.append(BedLine("\t".join(str(x) for x in d)))

    newbed.print_to_file(opts.outfile, sorted=True)


def density(args):
    """
    %prog density bedfile ref.fasta

    Calculates density of features per seqid.
    """
    p = OptionParser(density.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, fastafile = args
    bed = Bed(bedfile)
    sizes = Sizes(fastafile).mapping
    header = "seqid features size density_per_Mb".split()
    print("\t".join(header))
    for seqid, bb in bed.sub_beds():
        nfeats = len(bb)
        size = sizes[seqid]
        ds = nfeats * 1e6 / size
        print("\t".join(str(x) for x in (seqid, nfeats, size, "{0:.1f}".format(ds))))


def clr(args):
    """
    %prog clr [bamfile|bedpefile] ref.fasta

    Use mates from BEDPE to extract ranges where the ref is covered by mates.
    This is useful in detection of chimeric contigs.
    """
    p = OptionParser(clr.__doc__)
    p.set_bedpe()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedpe, ref = args
    if bedpe.endswith(".bam"):
        bedpefile = bedpe.replace(".bam", ".bedpe")
        if need_update(bedpe, bedpefile):
            cmd = "bamToBed -bedpe -i {0}".format(bedpe)
            sh(cmd, outfile=bedpefile)
        bedpe = bedpefile

    filtered = bedpe + ".filtered"
    if need_update(bedpe, filtered):
        filter_bedpe(
            bedpe, filtered, ref, rc=opts.rc, minlen=opts.minlen, maxlen=opts.maxlen
        )

    rmdup = filtered + ".sorted.rmdup"
    if need_update(filtered, rmdup):
        rmdup_bedpe(filtered, rmdup, dupwiggle=opts.dup)

    converted = rmdup + ".converted"
    if need_update(rmdup, converted):
        fp = open(rmdup)
        fw = open(converted, "w")
        for row in fp:
            r = BedpeLine(row)
            print(r.bedline, file=fw)
        fw.close()

    merged = converted + ".merge.bed"
    if need_update(converted, merged):
        mergeBed(converted)


def sfa_to_fq(sfa, qvchar):
    fq = sfa.rsplit(".", 1)[0] + ".fq"
    fp = must_open(sfa)
    fw = must_open(fq, "w")
    total = 0
    for row in fp:
        total += 1
        name, seq = row.split()
        qual = len(seq) * qvchar
        print("\n".join(("@" + name, seq, "+", qual)), file=fw)
    logging.debug("A total of {0} sequences written to `{1}`.".format(total, fq))
    return fq


def filter_bedpe(bedpe, filtered, ref, rc=False, rlen=None, minlen=2000, maxlen=8000):
    tag = " after RC" if rc else ""
    logging.debug(
        "Filter criteria: innie{0}, {1} <= insertsize <= {2}".format(
            tag, minlen, maxlen
        )
    )
    sizes = Sizes(ref).mapping
    fp = must_open(bedpe)
    fw = must_open(filtered, "w")
    retained = total = 0
    for row in fp:
        b = BedpeLine(row)
        total += 1
        if rc:
            b.rc()
        if not b.is_innie:
            continue
        b.score = b.outerdist
        if not minlen <= b.score <= maxlen:
            continue
        retained += 1
        if rlen:
            b.extend(rlen, sizes[b.seqid1])
        print(b, file=fw)
    logging.debug(
        "A total of {0} mates written to `{1}`.".format(
            percentage(retained, total), filtered
        )
    )
    fw.close()


def rmdup_bedpe(filtered, rmdup, dupwiggle=10):
    sortedfiltered = filtered + ".sorted"
    if need_update(filtered, sortedfiltered):
        sh("sort -k1,1 -k2,2n -i {0} -o {1}".format(filtered, sortedfiltered))

    logging.debug("Rmdup criteria: wiggle <= {0}".format(dupwiggle))
    fp = must_open(sortedfiltered)
    fw = must_open(rmdup, "w")
    data = [BedpeLine(x) for x in fp]
    retained = total = 0
    for seqid, ss in groupby(data, key=lambda x: x.seqid1):
        ss = list(ss)
        for i, a in enumerate(ss):
            if a.isdup:
                continue
            for b in ss[i + 1 :]:
                if b.start1 > a.start1 + dupwiggle:
                    break
                if b.isdup:
                    continue
                if (
                    a.seqid2 == b.seqid2
                    and a.start2 - dupwiggle <= b.start2 <= a.start2 + dupwiggle
                ):
                    b.isdup = True
        for a in ss:
            total += 1
            if a.isdup:
                continue
            retained += 1
            print(a, file=fw)
    logging.debug(
        "A total of {0} mates written to `{1}`.".format(
            percentage(retained, total), rmdup
        )
    )
    fw.close()


def alignextend(args):
    """
    %prog alignextend bedpefile ref.fasta

    Similar idea to alignextend, using mates from BEDPE and FASTA ref. See AMOS
    script here:

    https://github.com/nathanhaigh/amos/blob/master/src/Experimental/alignextend.pl
    """
    p = OptionParser(alignextend.__doc__)
    p.add_option("--len", default=100, type="int", help="Extend to this length")
    p.add_option(
        "--qv", default=31, type="int", help="Dummy qv score for extended bases"
    )
    p.add_option(
        "--bedonly",
        default=False,
        action="store_true",
        help="Only generate bed files, no FASTA",
    )
    p.set_bedpe()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedpe, ref = args
    qvchar = chr(opts.qv + 33)
    pf = bedpe.split(".")[0]

    filtered = bedpe + ".filtered"
    if need_update(bedpe, filtered):
        filter_bedpe(
            bedpe,
            filtered,
            ref,
            rc=opts.rc,
            minlen=opts.minlen,
            maxlen=opts.maxlen,
            rlen=opts.rlen,
        )

    rmdup = filtered + ".filtered.sorted.rmdup"
    if need_update(filtered, rmdup):
        rmdup_bedpe(filtered, rmdup, dupwiggle=opts.dup)

    if opts.bedonly:
        return

    bed1, bed2 = pf + ".1e.bed", pf + ".2e.bed"
    if need_update(rmdup, (bed1, bed2)):
        sh("cut -f1-3,7-9 {0}".format(rmdup), outfile=bed1)
        sh("cut -f4-6,7-8,10 {0}".format(rmdup), outfile=bed2)

    sfa1, sfa2 = pf + ".1e.sfa", pf + ".2e.sfa"
    if need_update((bed1, bed2, ref), (sfa1, sfa2)):
        for bed in (bed1, bed2):
            fastaFromBed(bed, ref, name=True, tab=True, stranded=True)

    fq1, fq2 = pf + ".1e.fq", pf + ".2e.fq"
    if need_update((sfa1, sfa2), (fq1, fq2)):
        for sfa in (sfa1, sfa2):
            sfa_to_fq(sfa, qvchar)


def seqids(args):
    """
    %prog seqids bedfile

    Print out all seqids on one line. Useful for graphics.karyotype.
    """
    p = OptionParser(seqids.__doc__)
    p.add_option("--maxn", default=100, type="int", help="Maximum number of seqids")
    p.add_option("--prefix", help="Seqids must start with")
    p.add_option("--exclude", default="random", help="Seqids should not contain")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    pf = opts.prefix
    exclude = opts.exclude
    bed = Bed(bedfile)
    s = bed.seqids
    if pf:
        s = [x for x in s if x.startswith(pf)]
    if exclude:
        s = [x for x in s if not exclude in x]
    s = s[: opts.maxn]
    print(",".join(s))


def juncs(args):
    """
    %prog junctions junctions1.bed [junctions2.bed ...]

    Given a TopHat junctions.bed file, trim the read overhang to get intron span

    If more than one junction bed file is provided, uniq the junctions and
    calculate cumulative (sum) junction support
    """
    from tempfile import mkstemp
    from pybedtools import BedTool

    p = OptionParser(juncs.__doc__)
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fh, trimbed = mkstemp(suffix=".bed")
    fw = must_open(trimbed, "w")
    for i, juncbed in enumerate(args):
        bed = Bed(juncbed, juncs=True)
        for b in bed:
            ovh = [int(x) for x in b.extra[-2].split(",")]
            b.start += ovh[0]
            b.end -= ovh[1]
            b.accn = "{0}-{1}".format(b.accn, i)
            b.extra = None
            print(b, file=fw)
    fw.close()

    if len(args) > 1:
        sh("sort -k1,1 -k2,2n {0} -o {0}".format(trimbed))

        tbed = BedTool(trimbed)
        grouptbed = tbed.groupby(g=[1, 2, 3, 6], c=5, ops=["sum"])

        cmd = """awk -F $'\t' 'BEGIN { OFS = FS } { ID = sprintf("mJUNC%07d", NR); print $1,$2,$3,ID,$5,$4; }'"""
        infile = grouptbed.fn
        sh(cmd, infile=infile, outfile=opts.outfile)
    else:
        sort([trimbed, "-o", opts.outfile])

    os.unlink(trimbed)


def random(args):
    """
    %prog random bedfile number_of_features

    Extract a random subset of features. Number of features can be an integer
    number, or a fractional number in which case a random fraction (for example
    0.1 = 10% of all features) will be extracted.
    """
    from random import sample
    from jcvi.formats.base import flexible_cast

    p = OptionParser(random.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, N = args
    assert is_number(N)

    b = Bed(bedfile)
    NN = flexible_cast(N)
    if NN < 1:
        NN = int(round(NN * len(b)))

    beds = sample(b, NN)
    new_bed = Bed()
    new_bed.extend(beds)

    outfile = bedfile.rsplit(".", 1)[0] + ".{0}.bed".format(N)
    new_bed.print_to_file(outfile)
    logging.debug("Write {0} features to `{1}`".format(NN, outfile))


def filter(args):
    """
    %prog filter bedfile

    Filter the bedfile to retain records between certain size range.
    """
    p = OptionParser(filter.__doc__)
    p.add_option("--minsize", default=0, type="int", help="Minimum feature length")
    p.add_option(
        "--maxsize", default=1000000000, type="int", help="Minimum feature length"
    )
    p.add_option(
        "--minaccn",
        type="int",
        help="Minimum value of accn, useful to filter based on coverage",
    )
    p.add_option("--minscore", type="int", help="Minimum score")
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    fp = must_open(bedfile)
    fw = must_open(opts.outfile, "w")
    minsize, maxsize = opts.minsize, opts.maxsize
    minaccn = opts.minaccn
    minscore = opts.minscore
    total = []
    keep = []
    for row in fp:
        try:
            b = BedLine(row)
        except IndexError:
            print(row.strip(), file=fw)
            continue
        span = b.span
        total.append(span)
        if not minsize <= span <= maxsize:
            continue
        if minaccn and int(b.accn) < minaccn:
            continue
        if minscore and int(b.score) < minscore:
            continue
        print(b, file=fw)
        keep.append(span)

    logging.debug("Stats: {0} features kept.".format(percentage(len(keep), len(total))))
    logging.debug("Stats: {0} bases kept.".format(percentage(sum(keep), sum(total))))


def make_bedgraph(bedfile, fastafile):
    sizesfile = Sizes(fastafile).filename
    pf = bedfile.rsplit(".", 1)[0]
    bedfile = sort([bedfile])
    bedgraph = pf + ".bedgraph"
    if need_update(bedfile, bedgraph):
        cmd = "genomeCoverageBed"
        cmd += " -i {0} -g {1} -bga".format(bedfile, sizesfile)
        sh(cmd, outfile=bedgraph)

    return bedgraph


def mergebydepth(args):
    """
    %prog mergebydepth reads.bed genome.fasta

    Similar to mergeBed, but only returns regions beyond certain depth.
    """
    p = OptionParser(mergebydepth.__doc__)
    p.add_option("--mindepth", default=3, type="int", help="Minimum depth required")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, fastafile = args
    mindepth = opts.mindepth
    bedgraph = make_bedgraph(bedfile)

    bedgraphfiltered = bedgraph + ".d{0}".format(mindepth)
    if need_update(bedgraph, bedgraphfiltered):
        filter(
            [
                bedgraph,
                "--minaccn={0}".format(mindepth),
                "--outfile={0}".format(bedgraphfiltered),
            ]
        )

    merged = bedgraphfiltered + ".merge.fasta"
    if need_update(bedgraphfiltered, merged):
        mergeBed(bedgraphfiltered, sorted=True)


def depth(args):
    """
    %prog depth reads.bed features.bed

    Calculate depth depth per feature using coverageBed.
    """
    p = OptionParser(depth.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    readsbed, featsbed = args
    fp = open(featsbed)
    nargs = len(fp.readline().split("\t"))
    keepcols = ",".join(str(x) for x in range(1, nargs + 1))
    cmd = "coverageBed -a {0} -b {1} -d".format(readsbed, featsbed)
    cmd += " | groupBy -g {0} -c {1} -o mean".format(keepcols, nargs + 2)
    sh(cmd, outfile=opts.outfile)


def remove_isoforms(ids):
    """
    This is more or less a hack to remove the GMAP multiple mappings. Multiple
    GMAP mappings can be seen given the names .mrna1, .mrna2, etc.
    """
    key = lambda x: x.rsplit(".", 1)[0]
    iso_number = lambda x: get_number(x.split(".")[-1])
    ids = sorted(ids, key=key)
    newids = []
    for k, ii in groupby(ids, key=key):
        min_i = min(list(ii), key=iso_number)
        newids.append(min_i)
    return newids


def longest(args):
    """
    %prog longest bedfile fastafile

    Select longest feature within overlapping piles.
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(longest.__doc__)
    p.add_option("--maxsize", default=20000, type="int", help="Limit max size")
    p.add_option("--minsize", default=60, type="int", help="Limit min size")
    p.add_option(
        "--precedence", default="Medtr", help="Accessions with prefix take precedence"
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, fastafile = args
    maxsize = opts.maxsize
    minsize = opts.minsize
    prec = opts.precedence
    mergedbed = mergeBed(bedfile, nms=True)
    sizes = Sizes(fastafile).mapping
    bed = Bed(mergedbed)

    pf = bedfile.rsplit(".", 1)[0]
    ids = set()
    for b in bed:
        accns = b.accn.split(";")
        prec_accns = [x for x in accns if x.startswith(prec)]
        if prec_accns:
            accns = prec_accns
        accn_sizes = [(sizes.get(x, 0), x) for x in accns]
        accn_sizes = [(size, x) for size, x in accn_sizes if size < maxsize]
        if not accn_sizes:
            continue
        max_size, max_accn = max(accn_sizes)
        if max_size < minsize:
            continue
        ids.add(max_accn)

    newids = remove_isoforms(ids)
    logging.debug("Remove isoforms: before={0} after={1}".format(len(ids), len(newids)))

    longestidsfile = pf + ".longest.ids"
    fw = open(longestidsfile, "w")
    print("\n".join(newids), file=fw)
    fw.close()
    logging.debug(
        "A total of {0} records written to `{1}`.".format(len(newids), longestidsfile)
    )

    longestbedfile = pf + ".longest.bed"
    some(
        [
            bedfile,
            longestidsfile,
            "--outfile={0}".format(longestbedfile),
            "--no_strip_names",
        ]
    )


def merge(args):
    """
    %prog merge bedfiles > newbedfile

    Concatenate bed files together. Performing seqid and name changes to avoid
    conflicts in the new bed file.
    """
    p = OptionParser(merge.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    bedfiles = args
    fw = must_open(opts.outfile, "w")
    for bedfile in bedfiles:
        bed = Bed(bedfile)
        pf = op.basename(bedfile).split(".")[0]
        for b in bed:
            b.seqid = "_".join((pf, b.seqid))
            print(b, file=fw)


def fix(args):
    """
    %prog fix bedfile > newbedfile

    Fix non-standard bed files. One typical problem is start > end.
    """
    p = OptionParser(fix.__doc__)
    p.add_option("--minspan", default=0, type="int", help="Enforce minimum span")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    minspan = opts.minspan
    fp = open(bedfile)
    fw = must_open(opts.outfile, "w")
    nfixed = nfiltered = ntotal = 0
    for row in fp:
        atoms = row.strip().split("\t")
        assert len(atoms) >= 3, "Must be at least 3 columns"
        seqid, start, end = atoms[:3]
        start, end = int(start), int(end)
        orientation = "+"
        if start > end:
            start, end = end, start
            orientation = "-"
            nfixed += 1

        atoms[1:3] = [str(start), str(end)]
        if len(atoms) > 6:
            atoms[6] = orientation
        line = "\t".join(atoms)
        b = BedLine(line)

        if b.span >= minspan:
            print(b, file=fw)
            nfiltered += 1

        ntotal += 1

    if nfixed:
        logging.debug("Total fixed: {0}".format(percentage(nfixed, ntotal)))
    if nfiltered:
        logging.debug("Total filtered: {0}".format(percentage(nfiltered, ntotal)))


def some(args):
    """
    %prog some bedfile idsfile > newbedfile

    Retrieve a subset of bed features given a list of ids.
    """
    from jcvi.formats.base import SetFile
    from jcvi.utils.cbook import gene_name

    p = OptionParser(some.__doc__)
    p.add_option(
        "-v",
        dest="inverse",
        default=False,
        action="store_true",
        help="Get the inverse, like grep -v",
    )
    p.set_outfile()
    p.set_stripnames()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, idsfile = args
    inverse = opts.inverse
    ostrip = opts.strip_names
    fw = must_open(opts.outfile, "w")

    ids = SetFile(idsfile)
    if ostrip:
        ids = set(gene_name(x) for x in ids)
    bed = Bed(bedfile)
    ntotal = nkeep = 0
    for b in bed:
        ntotal += 1
        keep = b.accn in ids
        if inverse:
            keep = not keep

        if keep:
            nkeep += 1
            print(b, file=fw)

    fw.close()
    logging.debug("Stats: {0} features kept.".format(percentage(nkeep, ntotal)))


def uniq(args):
    """
    %prog uniq bedfile

    Remove overlapping features with higher scores.
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(uniq.__doc__)
    p.add_option("--sizes", help="Use sequence length as score")
    p.add_option("--mode", default="span", choices=("span", "score"), help="Pile mode")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    uniqbedfile = bedfile.split(".")[0] + ".uniq.bed"
    bed = Bed(bedfile)

    if opts.sizes:
        sizes = Sizes(opts.sizes).mapping
        ranges = [
            Range(x.seqid, x.start, x.end, sizes[x.accn], i) for i, x in enumerate(bed)
        ]
    else:
        if opts.mode == "span":
            ranges = [
                Range(x.seqid, x.start, x.end, x.end - x.start + 1, i)
                for i, x in enumerate(bed)
            ]
        else:
            ranges = [
                Range(x.seqid, x.start, x.end, float(x.score), i)
                for i, x in enumerate(bed)
            ]

    selected, score = range_chain(ranges)
    selected = [x.id for x in selected]
    selected_ids = set(selected)
    selected = [bed[x] for x in selected]
    notselected = [x for i, x in enumerate(bed) if i not in selected_ids]

    newbed = Bed()
    newbed.extend(selected)
    newbed.print_to_file(uniqbedfile, sorted=True)

    if notselected:
        leftoverfile = bedfile.split(".")[0] + ".leftover.bed"
        leftoverbed = Bed()
        leftoverbed.extend(notselected)
        leftoverbed.print_to_file(leftoverfile, sorted=True)

    logging.debug("Imported: {0}, Exported: {1}".format(len(bed), len(newbed)))

    return uniqbedfile


def subtractbins(binfile1, binfile2):
    from jcvi.graphics.landscape import BinFile

    abin = BinFile(binfile1)
    bbin = BinFile(binfile2)

    assert len(abin) == len(bbin)

    fw = open(binfile1, "w")

    for a, b in zip(abin, bbin):
        assert a.chr == b.chr
        assert a.binlen == b.binlen

        a.subtract(b)
        print(a, file=fw)

    fw.close()

    return binfile1


def bins(args):
    """
    %prog bins bedfile fastafile

    Bin bed lengths into each consecutive window. Use --subtract to remove bases
    from window, e.g. --subtract gaps.bed ignores the gap sequences.
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(bins.__doc__)
    p.add_option("--binsize", default=100000, type="int", help="Size of the bins")
    p.add_option("--subtract", help="Subtract bases from window")
    p.add_option(
        "--mode",
        default="span",
        choices=("span", "count", "score"),
        help="Accumulate feature based on",
    )
    p.add_option(
        "--nomerge", default=False, action="store_true", help="Do not merge features"
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, fastafile = args
    subtract = opts.subtract
    mode = opts.mode
    assert op.exists(bedfile), "File `{0}` not found".format(bedfile)

    binsize = opts.binsize
    binfile = bedfile + ".{0}".format(binsize)
    binfile += ".{0}.bins".format(mode)

    if not need_update(bedfile, binfile):
        return binfile

    sz = Sizes(fastafile)
    sizesfile = sz.filename
    sizes = sz.mapping
    fw = open(binfile, "w")
    scores = "median" if mode == "score" else None
    if not opts.nomerge:
        bedfile = mergeBed(bedfile, nms=True, scores=scores)
    if subtract:
        subtractmerge = mergeBed(subtract)
        subtract_complement = complementBed(subtractmerge, sizesfile)
        bedfile = intersectBed(bedfile, subtract_complement)

    bedfile = sort([bedfile, "-i"])

    bed = Bed(bedfile)
    sbdict = dict(bed.sub_beds())
    for chr, chr_len in sorted(sizes.items()):
        chr_len = sizes[chr]
        subbeds = sbdict.get(chr, [])
        nbins = chr_len / binsize
        last_bin = chr_len % binsize
        if last_bin:
            nbins += 1

        a = np.zeros(nbins)  # values
        b = np.zeros(nbins, dtype="int")  # bases
        c = np.zeros(nbins, dtype="int")  # count
        b[:-1] = binsize
        b[-1] = last_bin

        for bb in subbeds:

            start, end = bb.start, bb.end
            startbin = start / binsize
            endbin = end / binsize

            assert startbin <= endbin
            c[startbin : endbin + 1] += 1

            if mode == "score":
                a[startbin : endbin + 1] += float(bb.score)

            elif mode == "span":
                if startbin == endbin:
                    a[startbin] += end - start + 1

                if startbin < endbin:
                    firstsize = (startbin + 1) * binsize - start + 1
                    lastsize = end - endbin * binsize
                    a[startbin] += firstsize
                    if startbin + 1 < endbin:
                        a[startbin + 1 : endbin] += binsize
                    a[endbin] += lastsize

        if mode == "count":
            a = c

        for xa, xb in zip(a, b):
            print("\t".join(str(x) for x in (chr, xa, xb)), file=fw)

    fw.close()

    if subtract:
        subtractbinfile = bins([subtract, fastafile, "--binsize={0}".format(binsize)])
        binfile = subtractbins(binfile, subtractbinfile)

    return binfile


def pile(args):
    """
    %prog pile abedfile bbedfile > piles

    Call intersectBed on two bedfiles.
    """
    from jcvi.utils.grouper import Grouper

    p = OptionParser(pile.__doc__)
    p.add_option("--minOverlap", default=0, type="int", help="Minimum overlap required")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    abedfile, bbedfile = args
    iw = intersectBed_wao(abedfile, bbedfile, minOverlap=opts.minOverlap)
    groups = Grouper()
    for a, b in iw:
        groups.join(a.accn, b.accn)

    ngroups = 0
    for group in groups:
        if len(group) > 1:
            ngroups += 1
            print("|".join(group))

    logging.debug("A total of {0} piles (>= 2 members)".format(ngroups))


def index(args):
    """
    %prog index bedfile

    Compress and index bedfile using `tabix`. Use --fasta to give a FASTA file
    so that a bedgraph file can be generated and indexed.
    """
    p = OptionParser(index.__doc__)
    p.add_option("--fasta", help="Generate bedgraph and index")
    p.add_option("--query", help="Chromosome location")
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    fastafile = opts.fasta
    if fastafile:
        bedfile = make_bedgraph(bedfile, fastafile)

    bedfile = sort([bedfile])

    gzfile = bedfile + ".gz"
    if need_update(bedfile, gzfile):
        cmd = "bgzip {0}".format(bedfile)
        sh(cmd)

    tbifile = gzfile + ".tbi"
    if need_update(gzfile, tbifile):
        cmd = "tabix -p bed {0}".format(gzfile)
        sh(cmd)

    query = opts.query
    if not query:
        return

    cmd = "tabix {0} {1}".format(gzfile, query)
    sh(cmd, outfile=opts.outfile)


def fastaFromBed(bedfile, fastafile, name=False, tab=False, stranded=False):
    suffix = ".sfa" if tab else ".fasta"
    outfile = op.basename(bedfile).rsplit(".", 1)[0] + suffix
    cmd = "fastaFromBed -fi {0} -bed {1} -fo {2}".format(fastafile, bedfile, outfile)
    if name:
        cmd += " -name"
    if tab:
        cmd += " -tab"
    if stranded:
        cmd += " -s"

    if need_update([bedfile, fastafile], outfile):
        sh(cmd, outfile=outfile)

    return outfile


def mergeBed(bedfile, d=0, sorted=False, nms=False, s=False, scores=None, delim=";"):
    if not sorted:
        bedfile = sort([bedfile, "-i"])
    cmd = "mergeBed -i {0}".format(bedfile)
    if d:
        cmd += " -d {0}".format(d)
    if nms:
        nargs = len(open(bedfile).readline().split())
        if nargs <= 3:
            logging.debug("Only {0} columns detected... set nms=True".format(nargs))
        else:
            cmd += " -c 4 -o collapse"
    if s:
        cmd += " -s"
    if scores:
        valid_opts = (
            "sum",
            "min",
            "max",
            "mean",
            "median",
            "mode",
            "antimode",
            "collapse",
        )
        if not scores in valid_opts:
            scores = "mean"
        cmd += " -scores {0}".format(scores)

    if delim:
        cmd += ' -delim "{0}"'.format(delim)

    pf = bedfile.rsplit(".", 1)[0] if bedfile.endswith(".bed") else bedfile
    mergebedfile = op.basename(pf) + ".merge.bed"

    if need_update(bedfile, mergebedfile):
        sh(cmd, outfile=mergebedfile)
    return mergebedfile


def complementBed(bedfile, sizesfile):
    cmd = "complementBed"
    cmd += " -i {0} -g {1}".format(bedfile, sizesfile)
    complementbedfile = "complement_" + op.basename(bedfile)

    if need_update([bedfile, sizesfile], complementbedfile):
        sh(cmd, outfile=complementbedfile)
    return complementbedfile


def intersectBed(bedfile1, bedfile2):
    cmd = "intersectBed"
    cmd += " -a {0} -b {1}".format(bedfile1, bedfile2)
    suffix = ".intersect.bed"

    intersectbedfile = (
        ".".join(
            (op.basename(bedfile1).split(".")[0], op.basename(bedfile2).split(".")[0])
        )
        + suffix
    )

    if need_update([bedfile1, bedfile2], intersectbedfile):
        sh(cmd, outfile=intersectbedfile)
    return intersectbedfile


def query_to_range(query, sizes):
    # chr1:1-10000 => (chr1, 0, 10000)
    if ":" in query:
        a, bc = query.split(":", 1)
        b, c = [int(x) for x in bc.split("-", 1)]
        b -= 1
    else:
        a = query
        b, c = 0, sizes.mapping[a]

    return a, b, c


def evaluate(args):
    """
    %prog evaluate prediction.bed reality.bed fastafile

    Make a truth table like:
            True    False  --- Reality
    True    TP      FP
    False   FN      TN
     |----Prediction

    Sn = TP / (all true in reality) = TP / (TP + FN)
    Sp = TP / (all true in prediction) = TP / (TP + FP)
    Ac = (TP + TN) / (TP + FP + FN + TN)
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(evaluate.__doc__)
    p.add_option("--query", help="Chromosome location")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    prediction, reality, fastafile = args
    query = opts.query
    prediction = mergeBed(prediction)
    reality = mergeBed(reality)
    sizes = Sizes(fastafile)
    sizesfile = sizes.filename

    prediction_complement = complementBed(prediction, sizesfile)
    reality_complement = complementBed(reality, sizesfile)

    TPbed = intersectBed(prediction, reality)
    FPbed = intersectBed(prediction, reality_complement)
    FNbed = intersectBed(prediction_complement, reality)
    TNbed = intersectBed(prediction_complement, reality_complement)
    beds = (TPbed, FPbed, FNbed, TNbed)

    if query:
        subbeds = []
        rr = query_to_range(query, sizes)
        ce = 'echo "{0}"'.format("\t".join(str(x) for x in rr))
        for b in beds:
            subbed = ".".join((b, query))
            cmd = ce + " | intersectBed -a stdin -b {0}".format(b)
            sh(cmd, outfile=subbed)
            subbeds.append(subbed)
        beds = subbeds

    be = BedEvaluate(*beds)
    print(be, file=sys.stderr)

    if query:
        for b in subbeds:
            os.remove(b)

    return be


def intersectBed_wao(abedfile, bbedfile, minOverlap=0):
    abed = Bed(abedfile)
    bbed = Bed(bbedfile)
    print("`{0}` has {1} features.".format(abedfile, len(abed)), file=sys.stderr)
    print("`{0}` has {1} features.".format(bbedfile, len(bbed)), file=sys.stderr)

    cmd = "intersectBed -wao -a {0} -b {1}".format(abedfile, bbedfile)
    acols = abed[0].nargs
    bcols = bbed[0].nargs
    fp = popen(cmd)
    for row in fp:
        atoms = row.split()
        aline = "\t".join(atoms[:acols])
        bline = "\t".join(atoms[acols : acols + bcols])
        c = int(atoms[-1])
        if c < minOverlap:
            continue
        a = BedLine(aline)
        try:
            b = BedLine(bline)
        except AssertionError:
            b = None

        yield a, b


def refine(args):
    """
    %prog refine bedfile1 bedfile2 refinedbed

    Refine bed file using a second bed file. The final bed is keeping all the
    intervals in bedfile1, but refined by bedfile2 whenever they have
    intersection.
    """
    p = OptionParser(refine.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    abedfile, bbedfile, refinedbed = args
    fw = open(refinedbed, "w")
    intersected = refined = 0
    for a, b in intersectBed_wao(abedfile, bbedfile):
        if b is None:
            print(a, file=fw)
            continue

        intersected += 1
        aspan_before = a.span
        arange = (a.start, a.end)
        brange = (b.start, b.end)
        irange = range_intersect(arange, brange)
        a.start, a.end = irange
        aspan_after = a.span
        if aspan_before > aspan_after:
            refined += 1
        print(a, file=fw)

    fw.close()
    print("Total intersected: {0}".format(intersected), file=sys.stderr)
    print("Total refined: {0}".format(refined), file=sys.stderr)
    summary([abedfile])
    summary([refinedbed])


def distance(args):
    """
    %prog distance bedfile

    Calculate distance between bed features. The output file is a list of
    distances, which can be used to plot histogram, etc.
    """
    p = OptionParser(distance.__doc__)
    p.add_option(
        "--distmode",
        default="ss",
        choices=("ss", "ee"),
        help="Distance mode between paired reads. ss is outer distance, "
        "ee is inner distance",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    sortedbedfile = sort([bedfile])
    valid = total = 0
    fp = open(sortedbedfile)
    for a, b in pairwise(fp):
        a = BedLine(a)
        b = BedLine(b)
        ar = (a.seqid, a.start, a.end, "+")
        br = (b.seqid, b.start, b.end, "+")
        dist, oo = range_distance(ar, br, distmode=opts.distmode)
        total += 1
        if dist > 0:
            print(dist)
            valid += 1

    logging.debug("Total valid (> 0) distances: {0}.".format(percentage(valid, total)))


def sample(args):
    """
    %prog sample bedfile sizesfile

    Sample bed file and remove high-coverage regions.

    When option --targetsize is used, this program uses a differnent mode. It
    first calculates the current total bases from all ranges and then compare to
    targetsize, if more, then sample down as close to targetsize as possible.

    Selection via --raindrop has the effect of making coverage even. Selected
    reads have the property that their end points are not within a certain
    window from one another. One sweep goes from left to right, the other in
    the reverse direction.
    """
    import random
    from jcvi.assembly.coverage import Coverage

    p = OptionParser(sample.__doc__)
    p.add_option(
        "--raindrop",
        default=0,
        type="int",
        help="Raindrop selection, ignores all other options",
    )
    p.add_option("--max", default=10, type="int", help="Max depth allowed")
    p.add_option(
        "--targetsize", type="int", help="Sample bed file to get target base number"
    )
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, sizesfile = args
    pf = bedfile.rsplit(".", 1)[0]
    raindrop = opts.raindrop

    # Raindrop method
    if raindrop:
        bed = Bed(bedfile)
        forward = []
        for b in bed:
            if not forward or abs(b.start - forward[-1].start) >= raindrop:
                forward.append(b)

        reverse = []
        bed.sort(key=lambda x: -x.end)
        for b in bed:
            if not reverse or abs(b.end - reverse[-1].end) >= raindrop:
                reverse.append(b)

        for tag, L in zip(("forward", "reverse"), (forward, reverse)):
            logging.debug(
                "Selected {0} features in {1} direction, span: {2}".format(
                    len(L), tag, sum(x.span for x in L)
                )
            )

        selected = Bed()
        selected.extend(set(forward + reverse))
        selected.print_to_file(opts.outfile, sorted=True)
        return

    targetsize = opts.targetsize
    if targetsize:
        bed = Bed(bedfile)
        samplebed = pf + ".sample.bed"
        fw = open(samplebed, "w")
        nfeats = len(bed)
        nbases = bed.sum(unique=False)
        targetfeats = int(round(nfeats * targetsize / nbases))
        sub_bed = random.sample(bed, targetfeats)
        for b in sub_bed:
            print(b, file=fw)

        logging.debug("File written to `%s`.", samplebed)
        return

    c = Coverage(bedfile, sizesfile)
    coveragefile = c.filename
    samplecoveragefile = pf + ".sample.coverage"
    fw = open(samplecoveragefile, "w")
    fp = open(coveragefile)
    for row in fp:
        seqid, start, end, cov = row.split()
        cov = int(cov)
        if cov <= opts.max:
            fw.write(row)
    fw.close()

    samplebedfile = pf + ".sample.bed"
    cmd = "intersectBed -a {0} -b {1} -wa -u".format(bedfile, samplecoveragefile)
    sh(cmd, outfile=samplebedfile)
    logging.debug("Sampled bedfile written to `{0}`.".format(samplebedfile))


def bedpe(args):
    """
    %prog bedpe bedfile

    Convert to bedpe format. Use --span to write another bed file that contain
    the span of the read pairs.
    """
    from jcvi.assembly.coverage import bed_to_bedpe

    p = OptionParser(bedpe.__doc__)
    p.add_option(
        "--span", default=False, action="store_true", help="Write span bed file"
    )
    p.add_option(
        "--strand", default=False, action="store_true", help="Write the strand columns"
    )
    p.add_option("--mates", help="Check the library stats from .mates file")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    pf = bedfile.rsplit(".", 1)[0]
    bedpefile = pf + ".bedpe"
    bedspanfile = pf + ".spans.bed" if opts.span else None
    bed_to_bedpe(
        bedfile,
        bedpefile,
        pairsbedfile=bedspanfile,
        matesfile=opts.mates,
        strand=opts.strand,
    )
    return bedpefile, bedspanfile


def sizes(args):
    """
    %prog sizes bedfile

    Infer the sizes for each seqid. Useful before dot plots.
    """
    p = OptionParser(sizes.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    assert op.exists(bedfile)

    sizesfile = bedfile.rsplit(".", 1)[0] + ".sizes"

    fw = must_open(sizesfile, "w", checkexists=True, skipcheck=True)
    if fw:
        b = Bed(bedfile)
        for s, sbeds in b.sub_beds():
            print("{0}\t{1}".format(s, max(x.end for x in sbeds)), file=fw)
        logging.debug("Sizes file written to `{0}`.".format(sizesfile))

    return sizesfile


def analyze_dists(dists, cutoff=1000, alpha=0.1):
    """
    The dists can show bimodal distribution if they come from a mate-pair
    library. Assume bimodal distribution and then separate the two peaks. Based
    on the percentage in each peak, we can decide if it is indeed one peak or
    two peaks, and report the median respectively.
    """
    peak0 = [d for d in dists if d < cutoff]
    peak1 = [d for d in dists if d >= cutoff]
    c0, c1 = len(peak0), len(peak1)
    logging.debug("Component counts: {0} {1}".format(c0, c1))
    if c0 == 0 or c1 == 0 or float(c1) / len(dists) < alpha:
        logging.debug(
            "Single peak identified ({0} / {1} < {2})".format(c1, len(dists), alpha)
        )
        return np.median(dists)

    peak0_median = np.median(peak0)
    peak1_median = np.median(peak1)
    logging.debug(
        "Dual peaks identified: {0}bp ({1}), {2}bp ({3}) (selected)".format(
            int(peak0_median), c0, int(peak1_median), c1
        )
    )

    return peak1_median


def report_pairs(
    data,
    cutoff=0,
    mateorientation=None,
    pairsfile=None,
    insertsfile=None,
    rclip=1,
    ascii=False,
    bins=20,
    distmode="ss",
    mpcutoff=1000,
):
    """
    This subroutine is used by the pairs function in blast.py and cas.py.
    Reports number of fragments and pairs as well as linked pairs
    """
    allowed_mateorientations = ("++", "--", "+-", "-+")

    if mateorientation:
        assert mateorientation in allowed_mateorientations

    num_fragments, num_pairs = 0, 0

    all_dist = []
    linked_dist = []
    # +- (forward-backward) is `innie`, -+ (backward-forward) is `outie`
    orientations = defaultdict(int)

    # clip how many chars from end of the read name to get pair name
    key = (lambda x: x.accn[:-rclip]) if rclip else (lambda x: x.accn)
    data.sort(key=key)

    if pairsfile:
        pairsfw = open(pairsfile, "w")
    if insertsfile:
        insertsfw = open(insertsfile, "w")

    for pe, lines in groupby(data, key=key):
        lines = list(lines)
        if len(lines) != 2:
            num_fragments += len(lines)
            continue

        num_pairs += 1
        a, b = lines

        asubject, astart, astop = a.seqid, a.start, a.end
        bsubject, bstart, bstop = b.seqid, b.start, b.end

        aquery, bquery = a.accn, b.accn
        astrand, bstrand = a.strand, b.strand

        dist, orientation = range_distance(
            (asubject, astart, astop, astrand),
            (bsubject, bstart, bstop, bstrand),
            distmode=distmode,
        )

        if dist >= 0:
            all_dist.append((dist, orientation, aquery, bquery))

    # select only pairs with certain orientations - e.g. innies, outies, etc.
    if mateorientation:
        all_dist = [x for x in all_dist if x[1] == mateorientation]

    # try to infer cutoff as twice the median until convergence
    if cutoff <= 0:
        dists = np.array([x[0] for x in all_dist], dtype="int")
        p0 = analyze_dists(dists, cutoff=mpcutoff)
        cutoff = int(2 * p0)  # initial estimate
        cutoff = int(math.ceil(cutoff / bins)) * bins
        logging.debug(
            "Insert size cutoff set to {0}, ".format(cutoff)
            + "use '--cutoff' to override"
        )

    for dist, orientation, aquery, bquery in all_dist:
        if dist > cutoff:
            continue
        if cutoff > 2 * mpcutoff and dist < mpcutoff:
            continue

        linked_dist.append(dist)
        if pairsfile:
            print("{0}\t{1}\t{2}".format(aquery, bquery, dist), file=pairsfw)
        orientations[orientation] += 1

    print(
        "{0} fragments, {1} pairs ({2} total)".format(
            num_fragments, num_pairs, num_fragments + num_pairs * 2
        ),
        file=sys.stderr,
    )

    s = SummaryStats(linked_dist, dtype="int")
    num_links = s.size

    meandist, stdev = s.mean, s.sd
    p0, p1, p2 = s.median, s.p1, s.p2

    print(
        "%d pairs (%.1f%%) are linked (cutoff=%d)"
        % (num_links, num_links * 100.0 / num_pairs, cutoff),
        file=sys.stderr,
    )
    print(
        "mean distance between mates: {0} +/- {1}".format(meandist, stdev),
        file=sys.stderr,
    )
    print("median distance between mates: {0}".format(p0), file=sys.stderr)
    print("95% distance range: {0} - {1}".format(p1, p2), file=sys.stderr)
    print("\nOrientations:", file=sys.stderr)

    orientation_summary = []
    for orientation, count in sorted(orientations.items()):
        o = "{0}:{1}".format(orientation, percentage(count, num_links, mode=1))
        orientation_summary.append(o.split()[0])
        print(o, file=sys.stderr)

    if insertsfile:
        from jcvi.graphics.histogram import histogram

        print("\n".join(str(x) for x in linked_dist), file=insertsfw)
        insertsfw.close()
        prefix = insertsfile.rsplit(".", 1)[0]
        if prefix > 10:
            prefix = prefix.split("-")[0]
        osummary = " ".join(orientation_summary)
        title = "{0} ({1}; median:{2} bp)".format(prefix, osummary, p0)
        histogram(
            insertsfile,
            vmin=0,
            vmax=cutoff,
            bins=bins,
            xlabel="Insertsize",
            title=title,
            ascii=ascii,
        )
        if op.exists(insertsfile):
            os.remove(insertsfile)

    return s


def pairs(args):
    """
    See __doc__ for OptionParser.set_pairs().
    """
    p = OptionParser(pairs.__doc__)
    p.set_pairs()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args

    basename = bedfile.split(".")[0]
    insertsfile = ".".join((basename, "inserts"))
    bedfile = sort([bedfile, "--accn"])

    fp = open(bedfile)
    data = [BedLine(row) for i, row in enumerate(fp) if i < opts.nrows]

    ascii = not opts.pdf
    return (
        bedfile,
        report_pairs(
            data,
            opts.cutoff,
            opts.mateorientation,
            pairsfile=opts.pairsfile,
            insertsfile=insertsfile,
            rclip=opts.rclip,
            ascii=ascii,
            bins=opts.bins,
            distmode=opts.distmode,
        ),
    )


def summary(args):
    """
    %prog summary bedfile

    Sum the total lengths of the intervals.
    """
    p = OptionParser(summary.__doc__)
    p.add_option(
        "--sizes", default=False, action="store_true", help="Write .sizes file"
    )
    p.add_option(
        "--all",
        default=False,
        action="store_true",
        help="Write summary stats per seqid",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    bed = Bed(bedfile)
    bs = BedSummary(bed)
    if opts.sizes:
        sizesfile = bedfile + ".sizes"
        fw = open(sizesfile, "w")
        for span, accn in bs.mspans:
            print(span, file=fw)
        fw.close()
        logging.debug("Spans written to `{0}`.".format(sizesfile))
        return bs

    if not opts.all:
        bs.report()
        return bs

    for seqid, subbeds in bed.sub_beds():
        bs = BedSummary(subbeds)
        print("\t".join((seqid, str(bs))))


def sort(args):
    """
    %prog sort bedfile

    Sort bed file to have ascending order of seqid, then start. It uses the
    `sort` command.
    """
    p = OptionParser(sort.__doc__)
    p.add_option(
        "-i",
        "--inplace",
        dest="inplace",
        default=False,
        action="store_true",
        help="Sort bed file in place",
    )
    p.add_option(
        "-u",
        dest="unique",
        default=False,
        action="store_true",
        help="Uniqify the bed file",
    )
    p.add_option(
        "--accn",
        default=False,
        action="store_true",
        help="Sort based on the accessions",
    )
    p.set_outfile(outfile=None)
    p.set_tmpdir()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    inplace = opts.inplace

    if not inplace and ".sorted." in bedfile:
        return bedfile

    sortedbed = opts.outfile
    if inplace:
        sortedbed = bedfile
    elif opts.outfile is None:
        pf, sf = op.basename(bedfile).rsplit(".", 1)
        sortedbed = pf + ".sorted." + sf

    sortopt = (
        "-k1,1 -k2,2n -k3,3n -k4,4" if not opts.accn else "-k4,4 -k1,1 -k2,2n -k3,3n"
    )
    cmd = "sort"
    if opts.tmpdir:
        cmd += " -T {0}".format(opts.tmpdir)
    if opts.unique:
        cmd += " -u"
    cmd += " {0} {1} -o {2}".format(sortopt, bedfile, sortedbed)

    if inplace or need_update(bedfile, sortedbed):
        sh(cmd)

    return sortedbed


def mates(args):
    """
    %prog mates bedfile

    Generate the mates file by inferring from the names.
    """
    p = OptionParser(mates.__doc__)
    p.add_option(
        "--lib",
        default=False,
        action="store_true",
        help="Output library information along with pairs",
    )
    p.add_option(
        "--nointra",
        default=False,
        action="store_true",
        help="Remove mates that are intra-scaffold",
    )
    p.add_option(
        "--prefix",
        default=False,
        action="store_true",
        help="Only keep links between IDs with same prefix",
    )
    p.set_mates()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    rclip = opts.rclip

    key = (lambda x: x.accn[:-rclip]) if rclip else (lambda x: x.accn)
    bed = Bed(bedfile, key=key)

    pf = bedfile.rsplit(".", 1)[0]
    matesfile = pf + ".mates"
    lib = pf if opts.lib else None
    fw = open(matesfile, "w")
    if lib:
        bedfile, stats = pairs(
            [bedfile, "--rclip={0}".format(rclip), "--cutoff={0}".format(opts.cutoff)]
        )
        sv = int(2 * stats.sd)
        mindist = max(stats.mean - sv, 1)
        maxdist = stats.mean + sv
        print("\t".join(str(x) for x in ("library", pf, mindist, maxdist)), file=fw)

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
        print("\t".join(pair), file=fw)

        print(a, file=fwm)
        print(b, file=fwm)

    logging.debug(
        "Discard {0} frags and write {1} pairs to `{2}` and `{3}`.".format(
            num_fragments, num_pairs, matesfile, matesbedfile
        )
    )

    fw.close()
    fwm.close()

    return matesfile, matesbedfile


def flanking(args):
    """
    %prog flanking bedfile [options]

    Get up to n features (upstream or downstream or both) flanking a given position.
    """
    from numpy import array, argsort

    p = OptionParser(flanking.__doc__)
    p.add_option(
        "--chrom",
        default=None,
        type="string",
        help="chrom name of the position in query. Make sure it matches bedfile.",
    )
    p.add_option(
        "--coord", default=None, type="int", help="coordinate of the position in query."
    )
    p.add_option(
        "-n", default=10, type="int", help="number of flanking features to get"
    )
    p.add_option(
        "--side",
        default="both",
        choices=("upstream", "downstream", "both"),
        help="which side to get flanking features",
    )
    p.add_option(
        "--max_d", default=None, type="int", help="features <= max_d away from position"
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if any([len(args) != 1, opts.chrom is None, opts.coord is None]):
        sys.exit(not p.print_help())

    (bedfile,) = args
    position = (opts.chrom, opts.coord)
    n, side, maxd = opts.n, opts.side, opts.max_d

    chrombed = Bed(bedfile).sub_bed(position[0])

    if side == "upstream":
        data = [
            (abs(f.start - position[1]), f) for f in chrombed if f.start <= position[1]
        ]
    elif side == "downstream":
        data = [
            (abs(f.start - position[1]), f) for f in chrombed if f.start >= position[1]
        ]
    else:
        data = [(abs(f.start - position[1]), f) for f in chrombed]

    if maxd:
        data = [f for f in data if f[0] <= maxd]

    n += 1  # not counting self
    n = min(n, len(data))
    distances, subbed = zip(*data)
    distances = array(distances)
    idx = argsort(distances)[:n]
    flankingbed = [f for (i, f) in enumerate(subbed) if i in idx]

    fw = must_open(opts.outfile, "w")
    for atom in flankingbed:
        print(str(atom), file=fw)

    return (position, flankingbed)


if __name__ == "__main__":
    main()
