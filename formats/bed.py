"""
Classes to handle the .bed files
"""

import os
import os.path as op
import sys
import math
import logging
import numpy as np

from collections import defaultdict
from itertools import groupby

from jcvi.formats.base import LineFile, must_open, is_number, get_number
from jcvi.utils.iter import pairwise
from jcvi.utils.cbook import SummaryStats, thousands, percentage
from jcvi.utils.natsort import natsort_key
from jcvi.utils.range import Range, range_union, range_chain, \
            range_distance, range_intersect
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, \
            need_update, popen


class BedLine(object):
    # the Bed format supports more columns. we only need
    # the first 4, but keep the information in 'extra'.
    __slots__ = ("seqid", "start", "end", "accn",
                 "extra", "score", "strand", "args", "nargs")

    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.nargs = nargs = len(args)
        self.seqid = args[0]
        self.start = int(args[1]) + 1
        self.end = int(args[2])
        assert self.start <= self.end, \
                "start={0} end={1}".format(self.start, self.end)
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

    def __getitem__(self, key):
        return getattr(self, key)

    @property
    def span(self):
        return self.end - self.start + 1

    @property
    def range(self):
        strand = self.strand or '+'
        return (self.seqid, self.start, self.end, strand)

    def gffline(self, type='match', source='default'):
        score = "." if not self.score or \
                (self.score and not is_number(self.score)) \
                else self.score
        strand = "." if not self.strand else self.strand
        row = "\t".join((self.seqid, source, type,
            str(self.start), str(self.end), score,
            strand, '.', 'ID=' + self.accn))
        return row


class Bed(LineFile):

    def __init__(self, filename=None, key=None, sorted=True, juncs=False):
        super(Bed, self).__init__(filename)

        # the sorting key provides some flexibility in ordering the features
        # for example, user might not like the lexico-order of seqid
        self.nullkey = lambda x: (natsort_key(x.seqid), x.start, x.accn)
        self.key = key or self.nullkey

        if not filename:
            return

        for line in must_open(filename):
            if line[0] == "#" or (juncs and line.startswith('track name')):
                continue
            self.append(BedLine(line))

        if sorted:
            self.sort(key=self.key)

    def print_to_file(self, filename="stdout", sorted=False):
        if sorted:
            self.sort(key=self.key)

        fw = must_open(filename, "w")
        for b in self:
            if b.start < 1:
                logging.error("Start < 1. Reset start for `{0}`.".format(b.accn))
                b.start = 1
            print >> fw, b
        fw.close()

    def sum(self, seqid=None, unique=True):
        return bed_sum(self, seqid=seqid, unique=unique)

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
    def simple_bed(self):
        return [(b.seqid, i) for (i, b) in enumerate(self)]

    @property
    def links(self):
        r = []
        for s, sb in self.sub_beds():
            for a, b in pairwise(sb):
                r.append(((a.accn, a.strand), (b.accn, b.strand)))
        return r

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


class BedEvaluate (object):

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

        msg += "\nSensitivity [TP / (TP + FN)]: {0:.1f} %\n".\
                format(self.sensitivity * 100)
        msg += "Specificity [TP / (TP + FP)]: {0:.1f} %\n".\
                format(self.specificity * 100)
        msg += "Accuracy [(TP + TN) / (TP + FP + FN + TN)]: {0:.1f} %".\
                format(self.accuracy * 100)
        return msg

    @property
    def sensitivity(self):
        if self.TP + self.FN == 0:
            return 0
        return self.TP * 1. / (self.TP + self.FN)

    @property
    def specificity(self):
        if self.TP + self.FP == 0:
            return 0
        return self.TP * 1. / (self.TP + self.FP)

    @property
    def accuracy(self):
        if self.TP + self.FP + self.FN + self.TN == 0:
            return 0
        return (self.TP + self.TN) * 1. / \
               (self.TP + self.FP + self.FN + self.TN)

    @property
    def score(self):
        return "|".join(("{0:.3f}".format(x) for x in \
                    (self.sensitivity, self.specificity, self.accuracy)))


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
        self.coverage = self.total_bases * 1. / self.unique_bases

    def report(self):
        print >> sys.stderr, "Total seqids: {0}".format(self.nseqids)
        print >> sys.stderr, "Total ranges: {0}".format(self.nfeats)
        print >> sys.stderr, "Total unique bases: {0} bp".format(thousands(self.unique_bases))
        print >> sys.stderr, "Total bases: {0} bp".format(thousands(self.total_bases))
        print >> sys.stderr, "Estimated coverage: {0:.1f}x".format(self.coverage)
        print >> sys.stderr, self.stats
        maxspan, maxaccn = max(self.mspans)
        minspan, minaccn = min(self.mspans)
        print >> sys.stderr, "Longest: {0} ({1})".format(maxaccn, maxspan)
        print >> sys.stderr, "Shortest: {0} ({1})".format(minaccn, minspan)

    def __str__(self):
        return "\t".join(str(x) for x in (self.nfeats, self.unique_bases))


def bed_sum(beds, seqid=None, unique=True):
    if seqid:
        ranges = [(x.seqid, x.start, x.end) for x in beds \
                    if x.seqid == seqid]
    else:
        ranges = [(x.seqid, x.start, x.end) for x in beds]

    unique_sum = range_union(ranges)
    raw_sum = sum(x.span for x in beds)
    return unique_sum if unique else raw_sum


def main():

    actions = (
        ('depth', 'calculate average depth per feature using coverageBed'),
        ('sort', 'sort bed file'),
        ('merge', 'merge bed files'),
        ('index', 'index bed file using tabix'),
        ('bins', 'bin bed lengths into each window'),
        ('summary', 'summarize the lengths of the intervals'),
        ('evaluate', 'make truth table and calculate sensitivity and specificity'),
        ('pile', 'find the ids that intersect'),
        ('pairs', 'estimate insert size between paired reads from bedfile'),
        ('mates', 'print paired reads from bedfile'),
        ('sizes', 'infer the sizes for each seqid'),
        ('uniq', 'remove overlapping features with higher scores'),
        ('longest', 'select longest feature within overlapping piles'),
        ('bedpe', 'convert to bedpe format'),
        ('distance', 'calculate distance between bed features'),
        ('sample', 'sample bed file and remove high-coverage regions'),
        ('refine', 'refine bed file using a second bed file'),
        ('flanking', 'get n flanking features for a given position'),
        ('some', 'get a subset of bed features given a list'),
        ('fix', 'fix non-standard bed files'),
        ('filter', 'filter the bedfile to retain records between size range'),
        ('random', 'extract a random subset of features'),
        ('juncs', 'trim junctions.bed overhang to get intron, merge multiple beds'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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

    fh, trimbed = mkstemp(suffix = ".bed")
    fw = must_open(trimbed, "w")
    for i, juncbed in enumerate(args):
        bed = Bed(juncbed, juncs=True)
        for b in bed:
            ovh = [int(x) for x in b.extra[-2].split(",")]
            b.start += ovh[0]
            b.end -= ovh[1]
            b.accn = "{0}-{1}".format(b.accn, i)
            b.extra = None
            print >> fw, b
    fw.close()

    if len(args) > 1:
        sh("sort -k1,1 -k2,2n {0} -o {0}".format(trimbed))

        tbed = BedTool(trimbed)
        grouptbed = tbed.groupby(g=[1,2,3,6], c=5, ops=['sum'])

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
    p.add_option("--minsize", default=0, type="int",
                 help="Minimum feature length")
    p.add_option("--maxsize", default=1000000000, type="int",
                 help="Minimum feature length")
    p.add_option("--minaccn", type="int",
                 help="Minimum value of accn, useful to filter based on coverage")
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    fp = must_open(bedfile)
    fw = must_open(opts.outfile, "w")
    minsize, maxsize = opts.minsize, opts.maxsize
    minaccn = opts.minaccn
    total = []
    keep = []
    for row in fp:
        b = BedLine(row)
        span = b.span
        total.append(span)
        if not minsize <= span <= maxsize:
            continue
        if minaccn and int(b.accn) < minaccn:
            continue
        print >> fw, b
        keep.append(span)

    logging.debug("Stats: {0} features kept.".\
                    format(percentage(len(keep), len(total))))
    logging.debug("Stats: {0} bases kept.".\
                    format(percentage(sum(keep), sum(total))))


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
    p.add_option("--maxsize", default=20000, type="int",
                 help="Limit max size")
    p.add_option("--minsize", default=60, type="int",
                 help="Limit min size")
    p.add_option("--precedence", default="Medtr",
                 help="Accessions with prefix take precedence")
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
    logging.debug("Remove isoforms: before={0} after={1}".\
                    format(len(ids), len(newids)))

    longestidsfile = pf + ".longest.ids"
    fw = open(longestidsfile, "w")
    print >> fw, "\n".join(newids)
    fw.close()
    logging.debug("A total of {0} records written to `{1}`.".\
                    format(len(newids), longestidsfile))

    longestbedfile = pf + ".longest.bed"
    some([bedfile, longestidsfile, "--outfile={0}".format(longestbedfile),
            "--no_strip_names"])


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
            print >> fw, b


def fix(args):
    """
    %prog fix bedfile > newbedfile

    Fix non-standard bed files. One typical problem is start > end.
    """
    p = OptionParser(fix.__doc__)
    p.add_option("--minspan", default=0, type="int",
                 help="Enforce minimum span [default: %default]")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    minspan = opts.minspan
    fp = open(bedfile)
    fw = must_open(opts.outfile, "w")
    nfixed = nfiltered = ntotal = 0
    for row in fp:
        atoms = row.strip().split("\t")
        assert len(atoms) >= 3, "Must be at least 3 columns"
        seqid, start, end = atoms[:3]
        start, end = int(start), int(end)
        orientation = '+'
        if start > end:
            start, end = end, start
            orientation = '-'
            nfixed += 1

        atoms[1:3] = [str(start), str(end)]
        if len(atoms) > 6:
            atoms[6] = orientation
        line = "\t".join(atoms)
        b = BedLine(line)

        if b.span >= minspan:
            print >> fw, b
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
    p.add_option("-v", dest="inverse", default=False, action="store_true",
                 help="Get the inverse, like grep -v [default: %default]")
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
            print >> fw, b

    fw.close()
    logging.debug("Stats: {0} features kept.".\
                    format(percentage(nkeep, ntotal)))


def uniq(args):
    """
    %prog uniq bedfile

    Remove overlapping features with higher scores.
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(uniq.__doc__)
    p.add_option("--sizes", help="Use sequence length as score")
    p.add_option("--mode", default="span", choices=("span", "score"),
                 help="Pile mode")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    uniqbedfile = bedfile.split(".")[0] + ".uniq.bed"
    bed = Bed(bedfile)

    if opts.sizes:
        sizes = Sizes(opts.sizes).mapping
        ranges = [Range(x.seqid, x.start, x.end, sizes[x.accn], i) \
                    for i, x in enumerate(bed)]
    else:
        if opts.mode == "span":
            ranges = [Range(x.seqid, x.start, x.end, x.end - x.start + 1, i) \
                        for i, x in enumerate(bed)]
        else:
            ranges = [Range(x.seqid, x.start, x.end, float(x.score), i) \
                        for i, x in enumerate(bed)]

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
        print >> fw, a

    fw.close()

    return binfile1


def bins(args):
    """
    %prog bins bedfile fastafile

    Bin bed lengths into each consecutive window. Use --subtract to remove bases
    from window, e.g. --subtract gaps.bed ignores the gap sequences.
    """
    import numpy as np

    from jcvi.formats.sizes import Sizes

    p = OptionParser(bins.__doc__)
    p.add_option("--binsize", default=100000, type="int",
                 help="Size of the bins [default: %default]")
    p.add_option("--subtract",
                 help="Subtract bases from window [default: %default]")
    p.add_option("--mode", default="span", choices=("span", "count", "score"),
                 help="Accumulate feature based on [default: %default]")
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
            c[startbin:endbin + 1] += 1

            if mode == "score":
                a[startbin:endbin + 1] += float(bb.score)

            elif mode == "span":
                if startbin == endbin:
                    a[startbin] += end - start + 1

                if startbin < endbin:
                    firstsize = (startbin + 1) * binsize - start + 1
                    lastsize = end - endbin * binsize
                    a[startbin] += firstsize
                    if startbin + 1 < endbin:
                        a[startbin + 1:endbin] += binsize
                    a[endbin] += lastsize

        if mode == "count":
            a = c

        for xa, xb in zip(a, b):
            print >> fw, "\t".join(str(x) for x in (chr, xa, xb))

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
    p.add_option("--minOverlap", default=0, type="int",
                 help="Minimum overlap required [default: %default]")
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
            print "|".join(group)

    logging.debug("A total of {0} piles (>= 2 members)".format(ngroups))


def index(args):
    """
    %prog index bedfile

    Compress frgscffile.sorted and index it using `tabix`.
    """
    p = OptionParser(index.__doc__)
    p.add_option("--query",
                 help="Chromosome location [default: %default]")
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    gzfile = bedfile + ".gz"

    if need_update(bedfile, gzfile):
        bedfile = sort([bedfile])
        cmd = "bgzip -c {0}".format(bedfile)
        sh(cmd, outfile=gzfile)

    tbifile = gzfile + ".tbi"

    if need_update(gzfile, tbifile):
        cmd = "tabix -p bed {0}".format(gzfile)
        sh(cmd)

    query = opts.query
    if not query:
        return

    cmd = "tabix {0} {1}".format(gzfile, query)
    sh(cmd, outfile=opts.outfile)


def fastaFromBed(bedfile, fastafile, name=False, stranded=False):
    outfile = op.basename(bedfile).rsplit(".", 1)[0] + ".fasta"
    cmd = "fastaFromBed -fi {0} -bed {1} -fo {2}".\
            format(fastafile, bedfile, outfile)
    if name:
        cmd += " -name"
    if stranded:
        cmd += " -s"

    if need_update([bedfile, fastafile], outfile):
        sh(cmd, outfile=outfile)

    return outfile


def mergeBed(bedfile, d=0, sorted=False, nms=False, s=False, scores=None):
    if not sorted:
        sort([bedfile, "-i"])
    cmd = "mergeBed -i {0}".format(bedfile)
    if d:
        cmd += " -d {0}".format(d)
    if nms:
        nargs = len(open(bedfile).readline().split())
        if nargs <= 3:
            logging.debug("Only {0} columns detected... set nms=True"\
                            .format(nargs))
        else:
            cmd += " -c 4 -o collapse"
    if s:
        cmd += " -s"
    if scores:
        valid_opts = ("sum", "min", "max", "mean", "median",
                "mode", "antimode", "collapse")
        if not scores in valid_opts:
            scores = "mean"
        cmd += " -scores {0}".format(scores)

    mergebedfile = op.basename(bedfile).rsplit(".", 1)[0] + ".merge.bed"

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

    intersectbedfile = ".".join((op.basename(bedfile1).split(".")[0],
            op.basename(bedfile2).split(".")[0])) + suffix

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
    p.add_option("--query",
                 help="Chromosome location [default: %default]")
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
    print >> sys.stderr, be

    if query:
        for b in subbeds:
            os.remove(b)

    return be


def intersectBed_wao(abedfile, bbedfile, minOverlap=0):
    abed = Bed(abedfile)
    bbed = Bed(bbedfile)
    print >> sys.stderr, "`{0}` has {1} features.".format(abedfile, len(abed))
    print >> sys.stderr, "`{0}` has {1} features.".format(bbedfile, len(bbed))

    cmd = "intersectBed -wao -a {0} -b {1}".format(abedfile, bbedfile)
    acols = abed[0].nargs
    bcols = bbed[0].nargs
    fp = popen(cmd)
    for row in fp:
        atoms = row.split()
        aline = "\t".join(atoms[:acols])
        bline = "\t".join(atoms[acols:acols + bcols])
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
            print >> fw, a
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
        print >> fw, a

    fw.close()
    print >> sys.stderr, "Total intersected: {0}".format(intersected)
    print >> sys.stderr, "Total refined: {0}".format(refined)
    summary([abedfile])
    summary([refinedbed])


def distance(args):
    """
    %prog distance bedfile

    Calculate distance between bed features. The output file is a list of
    distances, which can be used to plot histogram, etc.
    """
    from jcvi.utils.iter import pairwise

    p = OptionParser(distance.__doc__)
    p.add_option("--distmode", default="ss", choices=("ss", "ee"),
            help="Distance mode between paired reads. ss is outer distance, " \
                 "ee is inner distance [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
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
            print dist
            valid += 1

    logging.debug("Total valid (> 0) distances: {0}.".\
                  format(percentage(valid, total)))


def sample(args):
    """
    %prog sample bedfile sizesfile

    Sample bed file and remove high-coverage regions.

    When option --targetsize is used, this program uses a differnent mode. It
    first calculates the current total bases from all ranges and then compare to
    targetsize, if more, then sample down as close to targetsize as possible.
    """
    import random
    from jcvi.assembly.coverage import Coverage

    p = OptionParser(sample.__doc__)
    p.add_option("--max", default=10, type="int",
                 help="Max depth allowed [default: %default]")
    p.add_option("--targetsize", type="int",
                 help="Sample bed file to get target base number [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, sizesfile = args
    pf = bedfile.rsplit(".", 1)[0]

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
            print >> fw, b

        logging.debug("File written to `{0}`.".format(samplebed))
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
    p.add_option("--span", default=False, action="store_true",
                 help="Write span bed file [default: %default]")
    p.add_option("--mates", help="Check the library stats from .mates file")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    pf = bedfile.rsplit(".", 1)[0]
    bedpefile = pf + ".bedpe"
    bedspanfile = pf + ".spans.bed" if opts.span else None
    bed_to_bedpe(bedfile, bedpefile, \
                 pairsbedfile=bedspanfile, matesfile=opts.mates)
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


def analyze_dists(dists, cutoff=1000, alpha=.1):
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
        logging.debug("Single peak identified ({0} / {1} < {2})".\
                        format(c1, len(dists), alpha))
        return np.median(dists)

    peak0_median = np.median(peak0)
    peak1_median = np.median(peak1)
    logging.debug("Dual peaks identified: {0}bp ({1}), {2}bp ({3}) (selected)".\
                format(int(peak0_median), c0, int(peak1_median), c1))

    return peak1_median


def report_pairs(data, cutoff=0, mateorientation=None,
        pairsfile=None, insertsfile=None, rclip=1, ascii=False, bins=20,
        distmode="ss", mpcutoff=1000):
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

        dist, orientation = range_distance(\
                (asubject, astart, astop, astrand),
                (bsubject, bstart, bstop, bstrand),
                distmode=distmode)

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
        logging.debug("Insert size cutoff set to {0}, ".format(cutoff) +
            "use '--cutoff' to override")

    for dist, orientation, aquery, bquery in all_dist:
        if dist > cutoff:
            continue
        if cutoff > 2 * mpcutoff and dist < mpcutoff:
            continue

        linked_dist.append(dist)
        if pairsfile:
            print >> pairsfw, "{0}\t{1}\t{2}".format(aquery, bquery, dist)
        orientations[orientation] += 1

    print >>sys.stderr, "{0} fragments, {1} pairs ({2} total)".\
                format(num_fragments, num_pairs, num_fragments + num_pairs * 2)

    s = SummaryStats(linked_dist, dtype="int")
    num_links = s.size

    meandist, stdev = s.mean, s.sd
    p0, p1, p2 = s.median, s.p1, s.p2

    print >>sys.stderr, "%d pairs (%.1f%%) are linked (cutoff=%d)" % \
            (num_links, num_links * 100. / num_pairs, cutoff)
    print >>sys.stderr, "mean distance between mates: {0} +/- {1}".\
            format(meandist, stdev)
    print >>sys.stderr, "median distance between mates: {0}".format(p0)
    print >>sys.stderr, "95% distance range: {0} - {1}".format(p1, p2)
    print >>sys.stderr, "\nOrientations:"

    orientation_summary = []
    for orientation, count in sorted(orientations.items()):
        o = "{0}:{1}".format(orientation, \
                percentage(count, num_links, mode=1))
        orientation_summary.append(o.split()[0])
        print >>sys.stderr, o

    if insertsfile:
        from jcvi.graphics.histogram import histogram

        print >>insertsfw, "\n".join(str(x) for x in linked_dist)
        insertsfw.close()
        prefix = insertsfile.rsplit(".", 1)[0]
        if prefix > 10:
            prefix = prefix.split("-")[0]
        osummary = " ".join(orientation_summary)
        title="{0} ({1}; median:{2} bp)".format(prefix, osummary, p0)
        histogram(insertsfile, vmin=0, vmax=cutoff, bins=bins,
                xlabel="Insertsize", title=title, ascii=ascii)
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

    bedfile, = args

    basename = bedfile.split(".")[0]
    insertsfile = ".".join((basename, "inserts"))

    sortedbedfile = op.basename(bedfile).rsplit(".", 1)[0] + ".sorted.bed"
    if need_update(bedfile, sortedbedfile):
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


def summary(args):
    """
    %prog summary bedfile

    Sum the total lengths of the intervals.
    """
    p = OptionParser(summary.__doc__)
    p.add_option("--sizes", default=False, action="store_true",
                 help="Write .sizes file")
    p.add_option("--all", default=False, action="store_true",
                 help="Write summary stats per seqid")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    bed = Bed(bedfile)
    bs = BedSummary(bed)
    if opts.sizes:
        sizesfile = bedfile + ".sizes"
        fw = open(sizesfile, "w")
        for span, accn in bs.mspans:
            print >> fw, span
        fw.close()
        logging.debug("Spans written to `{0}`.".format(sizesfile))
        return bs

    if not opts.all:
        bs.report()
        return bs

    for seqid, subbeds in bed.sub_beds():
        bs = BedSummary(subbeds)
        print "\t".join((seqid, str(bs)))


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
    p.add_option("-u", dest="unique",
            default=False, action="store_true",
            help="Uniqify the bed file")
    p.add_option("--accn", default=False, action="store_true",
            help="Sort based on the accessions [default: %default]")
    p.set_outfile(outfile=None)
    p.set_tmpdir()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    inplace = opts.inplace

    sortedbed = opts.outfile
    if inplace:
        sortedbed = bedfile
    elif opts.outfile is None:
        sortedbed = op.basename(bedfile).rsplit(".", 1)[0] + ".sorted.bed"

    sortopt = "-k1,1 -k2,2n -k4,4" if not opts.accn else \
              "-k4,4 -k1,1 -k2,2n"
    cmd = "sort"
    if opts.tmpdir:
        cmd += " -T {0}".format(opts.tmpdir)
    if opts.unique:
        cmd += " -u"
    cmd += " {0} {1} -o {2}".format(sortopt, bedfile, sortedbed)
    sh(cmd)

    return sortedbed


def mates(args):
    """
    %prog mates bedfile

    Generate the mates file by inferring from the names.
    """
    p = OptionParser(mates.__doc__)
    p.add_option("--lib", default=False, action="store_true",
            help="Output library information along with pairs [default: %default]")
    p.add_option("--nointra", default=False, action="store_true",
            help="Remove mates that are intra-scaffold [default: %default]")
    p.add_option("--prefix", default=False, action="store_true",
            help="Only keep links between IDs with same prefix [default: %default]")
    p.set_mates()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    rclip = opts.rclip

    key = (lambda x: x.accn[:-rclip]) if rclip else (lambda x: x.accn)
    bed = Bed(bedfile, key=key)

    pf = bedfile.rsplit(".", 1)[0]
    matesfile = pf + ".mates"
    lib = pf if opts.lib else None
    fw = open(matesfile, "w")
    if lib:
        bedfile, stats = pairs([bedfile, \
                        "--rclip={0}".format(rclip),
                        "--cutoff={0}".format(opts.cutoff)])
        sv = int(2 * stats.sd)
        mindist = max(stats.mean - sv, 1)
        maxdist = stats.mean + sv
        print >> fw, "\t".join(str(x) for x in \
                ("library", pf, mindist, maxdist))

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


def flanking(args):
    """
    %prog flanking bedfile [options]

    Get up to n features (upstream or downstream or both) flanking a given position.
    """
    from numpy import array, argsort

    p = OptionParser(flanking.__doc__)
    p.add_option("--chrom", default=None, type="string",
            help="chrom name of the position in query. Make sure it matches bedfile.")
    p.add_option("--coord", default=None, type="int",
            help="coordinate of the position in query.")
    p.add_option("-n", default=10, type="int",
            help="number of flanking features to get [default: %default]")
    p.add_option("--side", default="both", choices=("upstream", "downstream", "both"),
            help="which side to get flanking features [default: %default]")
    p.add_option("--max_d", default=None, type="int",
            help="features <= max_d away from position [default: %default]")
    p.set_outfile()

    opts, args = p.parse_args(args)

    if any([len(args) != 1, opts.chrom is None, opts.coord is None]):
        sys.exit(not p.print_help())

    bedfile, = args
    position = (opts.chrom, opts.coord)
    n, side, maxd = opts.n, opts.side, opts.max_d

    chrombed = Bed(bedfile).sub_bed(position[0])

    if side == "upstream":
        data = [(abs(f.start-position[1]), f) for f in chrombed \
            if f.start <= position[1]]
    elif side == "downstream":
        data = [(abs(f.start-position[1]), f) for f in chrombed \
            if f.start >= position[1]]
    else:
        data = [(abs(f.start-position[1]), f) for f in chrombed]

    if maxd:
        data = [f for f in data if f[0]<=maxd]

    n += 1 # not counting self
    n = min(n, len(data))
    distances, subbed = zip(*data)
    distances = array(distances)
    idx = argsort(distances)[:n]
    flankingbed = [f for (i, f) in enumerate(subbed) if i in idx]

    fw = must_open(opts.outfile, "w")
    for atom in flankingbed:
        print >>fw, str(atom)

    return (position, flankingbed)


if __name__ == '__main__':
    main()
