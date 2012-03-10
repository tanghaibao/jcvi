"""
Classes to handle the .bed files
"""

import os
import os.path as op
import sys
import shutil
import logging

from itertools import groupby
from optparse import OptionParser

from jcvi.formats.base import LineFile, must_open
from jcvi.utils.cbook import depends, thousands
from jcvi.utils.range import range_union
from jcvi.apps.base import ActionDispatcher, debug, sh, \
        need_update, popen, set_outfile
debug()


class BedLine(object):
    # the Bed format supports more columns. we only need
    # the first 4, but keep the information in 'extra'.
    __slots__ = ("seqid", "start", "end", "accn",
                 "extra", "score", "strand", "nargs")

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
            self.extra = args[4:]
            self.score = self.extra[0]
        if nargs > 5:
            self.strand = self.extra[1]

    def __str__(self):
        args = [self.seqid, self.start - 1, self.end]
        if self.accn:
            args += [self.accn]
        if self.extra:
            args += self.extra

        s = "\t".join(str(x) for x in args)
        return s

    def __getitem__(self, key):
        return getattr(self, key)

    @property
    def span(self):
        return self.end - self.start + 1

    def reverse_complement(self, sizes):
        # this function is used in assembly.sopra
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

    def gffline(self, type='match', source='default'):
        score = "1000" if self.score == '.' else self.score
        row = "\t".join((self.seqid, source, type,
            str(self.start + 1), str(self.end), score,
            self.strand, '.', 'ID=' + self.accn))
        return row


class Bed(LineFile):

    def __init__(self, filename=None, key=None):
        super(Bed, self).__init__(filename)

        # the sorting key provides some flexibility in ordering the features
        # for example, user might not like the lexico-order of seqid
        self.nullkey = lambda x: (x.seqid, x.start, x.accn)
        self.key = key or self.nullkey

        if not filename:
            return

        for line in open(filename):
            if line[0] == "#":
                continue
            self.append(BedLine(line))

        self.sort(key=self.key)

    def print_to_file(self, filename):
        fw = must_open(filename, "w")
        for bedline in self:
            print >> fw, bedline
        fw.close()

    def sum(self, seqid=None, unique=True):
        if seqid:
            ranges = [(x.seqid, x.start, x.end) for x in self \
                        if x.seqid == seqid]
        else:
            ranges = [(x.seqid, x.start, x.end) for x in self]

        unique_sum = range_union(ranges)
        raw_sum = sum(x.span for x in self)
        return unique_sum if unique else raw_sum

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


def main():

    actions = (
        ('sort', 'sort bed file'),
        ('index', 'index bed file using tabix'),
        ('summary', 'summarize the lengths of the intervals'),
        ('evaluate', 'make truth table and calculate sensitivity and specificity'),
        ('pairs', 'estimate insert size between paired reads from bedfile'),
        ('mates', 'print paired reads from bedfile'),
        ('sizes', 'infer the sizes for each seqid'),
        ('bedpe', 'convert to bedpe format'),
        ('distance', 'calculate distance between bed features'),
        ('sample', 'sample bed file and remove high-coverage regions'),
        ('refine', 'refine bed file using a second bed file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def index(args):
    """
    %prog index bedfile

    Compress frgscffile.sorted and index it using `tabix`.
    """
    p = OptionParser(index.__doc__)
    p.add_option("--query",
                 help="Chromosome location [default: %default]")
    set_outfile(p)

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


def mergeBed(bedfile):
    cmd = "mergeBed -i {0}".format(bedfile)
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
    intersectbedfile = ".".join((op.basename(bedfile1).split(".")[0],
            op.basename(bedfile2).split(".")[0])) + ".intersect.bed"

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


def intersectBed_wao(abedfile, bbedfile):
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
    from jcvi.utils.range import range_intersect

    p = OptionParser(refine.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    abedfile, bbedfile, refinedbed = args
    fw = open(refinedbed, "w")
    intersected = refined = 0
    for a, b in intersected_wao(abedfile, bbedfile):
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


def distance(args):
    """
    %prog distance bedfile

    Calculate distance between bed features. The output file is a list of
    distances, which can be used to plot histogram, etc.
    """
    from jcvi.utils.cbook import percentage
    from jcvi.utils.iter import pairwise
    from jcvi.utils.range import range_distance

    p = OptionParser(distance.__doc__)
    p.add_option("--distmode", default="ss", choices=("ss", "ee"),
            help="distance mode between paired reads, ss is outer distance, " \
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
    import numpy as np

    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    bedfile, = args
    bed = Bed(bedfile)
    spans = np.array([x.span for x in bed])
    avg = int(np.average(spans))
    std = int(np.std(spans))
    print >> sys.stderr, "Total seqids: {0}".format(len(bed.seqids))
    print >> sys.stderr, "Total ranges: {0}".format(len(bed))

    total_bases = bed.sum(unique=False)
    unique_bases = bed.sum()

    print >> sys.stderr, "Total unique bases: {0} bp".format(thousands(unique_bases))
    print >> sys.stderr, "Total bases: {0} bp".format(thousands(total_bases))
    print >> sys.stderr, "Estimated coverage: {0:.1f}x".\
                        format(total_bases * 1. / unique_bases)

    print >> sys.stderr, "Average spans: {0}, stdev: {1}".format(avg, std)


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
    p.add_option("--lib", default=False, action="store_true",
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

    key = (lambda x: x.accn[:-rclip]) if rclip else (lambda x: x.accn)
    bed = Bed(bedfile, key=key)

    pf = bedfile.rsplit(".", 1)[0]
    matesfile = pf + ".mates"
    lib = pf if opts.lib else None
    fw = open(matesfile, "w")
    if lib:
        bedfile, (meandist, stdev, p0, p1, p2) = pairs([bedfile, \
                "--mateorientation=+-"])
        sv = int(2 * stdev)
        mindist = max(meandist - sv, 1)
        maxdist = meandist + sv
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


if __name__ == '__main__':
    main()
