#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Deals with K-mers and K-mer distribution from reads or genome
"""

import os.path as op
import sys
import logging
import numpy as np

from jcvi.graphics.base import plt, asciiplot, set_human_axis, savefig, \
            markup, panel_labels, normalize_axes, set_ticklabels_helvetica
from jcvi.formats.fasta import Fasta
from jcvi.formats.base import BaseFile, must_open, get_number
from jcvi.utils.cbook import thousands, percentage
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, \
            need_update, Popen, PIPE


KMERYL, KSOAP, KALLPATHS = range(3)


class KmerSpectrum (BaseFile):

    def __init__(self, histfile):
        self.load_data(histfile)

    def load_data(self, histfile):
        self.data = []
        self.totalKmers = 0
        self.hist = {}
        kformat = self.guess_format(histfile)
        kformats = ("Meryl", "Soap", "AllPaths")
        logging.debug("Guessed format: {0}".format(kformats[kformat]))

        fp = open(histfile)
        for rowno, row in enumerate(fp):
            if row[0] == '#':
                continue
            if kformat == KSOAP:
                K = rowno + 1
                counts = int(row.strip())
            else:  # meryl histogram
                K, counts = row.split()[:2]
                K, counts = int(K), int(counts)

            Kcounts = K * counts
            self.totalKmers += Kcounts
            self.hist[K] = Kcounts
            self.data.append((K, counts))

    def guess_format(self, histfile):
        # Guess the format of the Kmer histogram
        fp = open(histfile)
        for row in fp:
            if row.startswith("# 1:"):
                return KALLPATHS
            if len(row.split()) == 1:
                return KSOAP
        return KMERYL

    def get_xy(self, vmin=1, vmax=100):
        self.counts = sorted((a, b) for a, b in self.hist.items() \
                        if vmin <= a <= vmax)
        return zip(*self.counts)

    def analyze(self, ploidy=2, K=23, covmax=1000000):
        """
        Analyze Kmer spectrum, calculations derived from
        allpathslg/src/kmers/KmerSpectra.cc
        """
        from math import sqrt

        data = self.data
        kf_ceil = max(K for (K, c) in data)
        if kf_ceil > covmax:
            exceeds = sum(1 for (K, c) in data if K > covmax)
            logging.debug("A total of {0} distinct K-mers appear > "
                          "{1} times. Ignored ...".format(exceeds, covmax))
            kf_ceil = covmax

        nkf = kf_ceil + 1
        a = [0] * nkf
        for kf, c in data:
            if kf > kf_ceil:
                continue
            a[kf] = c

        ndk = a  # number of distinct kmers
        nk = [k * c for k, c in enumerate(a)]  # number of kmers
        cndk = [0] * nkf  # cumulative number of distinct kmers
        cnk = [0] * nkf  # cumulative number of kmers
        for kf in xrange(1, nkf):
            cndk[kf] = cndk[kf - 1] + .5 * (ndk[kf - 1] + ndk[kf])
            cnk [kf] = cnk [kf - 1] + .5 * (nk [kf - 1] + nk [kf])

        # Separate kmer spectrum in 5 regions based on the kf
        # 1        ... kf_min1    : bad kmers with low frequency
        # kf_min1  ... kf_min2    : good kmers CN = 1/2 (SNPs)
        # kf_min2  ... kf_min3    : good kmers CN = 1
        # kf_min3  ... kf_hi      : good kmers CN > 1 (repetitive)
        # kf_hi    ... inf        : bad kmers with high frequency

        # min1: find first minimum
        _kf_min1 = 10
        while (_kf_min1 - 1 >= 2 and nk[_kf_min1 - 1] < nk[_kf_min1]):
            _kf_min1 -= 1
        while (_kf_min1 <= kf_ceil and nk [_kf_min1 + 1] < nk[_kf_min1]):
            _kf_min1 += 1

        # max2: find absolute maximum mx2 above first minimum min1
        _kf_max2 = _kf_min1
        for kf in xrange(_kf_min1 + 1, int(0.8 * kf_ceil)):
            if nk[kf] > nk[_kf_max2]:
                _kf_max2 = kf

        # max2: resetting max2 for cases of very high polymorphism
        if ploidy == 2:
            ndk_half = ndk[_kf_max2 / 2]
            ndk_double = ndk[_kf_max2 * 2]
            if ndk_double > ndk_half:
                _kf_max2 *= 2

        # max1: SNPs local maximum max1 as half global maximum max2
        _kf_max1 = _kf_max2 / 2

        # min2: SNPs local minimum min2 between max1 and max2
        _kf_min2 = _kf_max1 * (2 * ndk[_kf_max1] + ndk[_kf_max2]) / \
                    (ndk[_kf_max1] + ndk[_kf_max2])

        # min1: refine between min1 and max2/2
        for kf in xrange(_kf_min1 + 1, _kf_max1):
            if (nk[kf] < nk[_kf_min1]):
                _kf_min1 = kf

        # min3: not a minimum, really. upper edge of main peak
        _kf_min3 = _kf_max2 * 3 / 2

        print >> sys.stderr, "kfs:", _kf_min1, _kf_max1, \
                        _kf_min2, _kf_max2, _kf_min3
        self.min1 = _kf_min1
        self.max1 = _kf_max1
        self.min2 = _kf_min2
        self.max2 = _kf_max2
        self.min3 = _kf_min3

        # Define maximum kf above which we neglect data
        _kf_hi = _kf_max2 * sqrt(4 * ndk[2 * _kf_max2] * _kf_max2) \
                    if 2 * _kf_max2 < len(ndk) else \
                 _kf_max2 * sqrt(4 * ndk[len(ndk) - 1] * _kf_max2)
        _kf_hi = int(_kf_hi)

        if _kf_hi > kf_ceil:
            _kf_hi = kf_ceil

        _nk_total = cnk[len(cnk) - 1]
        _nk_bad_low_kf = cnk[_kf_min1]
        _nk_good_uniq = cnk[_kf_min3] - cnk[_kf_min2]
        _nk_bad_high_kf = _nk_total - cnk[_kf_hi]
        _ndk_good_snp = cndk[_kf_min2] - cndk[_kf_min1]
        _ndk_good_uniq = cndk[_kf_min3] - cndk[_kf_min2]

        # kmer coverage C_k
        _kf_ave_uniq = _nk_good_uniq * 1. / _ndk_good_uniq
        _genome_size = (_nk_total - _nk_bad_low_kf - _nk_bad_high_kf) / \
                        _kf_ave_uniq
        _genome_size_unique = _ndk_good_uniq + _ndk_good_snp / 2
        _genome_size_repetitive = _genome_size - _genome_size_unique
        _coverage = _nk_total / _genome_size if _genome_size else 0

        # SNP rate estimation, assumes uniform distribution of SNPs over the
        # genome and accounts for the reduction in SNP kmer counts when
        # polymorphism is very high
        if ploidy == 2:
            _d_SNP = 1. / (1. - (1. - .5 * _ndk_good_snp / _genome_size) ** (1. / K)) \
                     if _ndk_good_snp > 0 else 1000000

        G = int(_genome_size)
        G1 = int(_genome_size_unique)
        GR = int(_genome_size_repetitive)
        coverage = int(_coverage)

        m = "Kmer (K={0}) Spectrum Analysis\n".format(K)
        m += "Genome size estimate = {0}\n".format(thousands(G))
        m += "Genome size estimate CN = 1 = {0} ({1})\n".format(thousands(G1),
                percentage(G1, G))
        m += "Genome size estimate CN > 1 = {0} ({1})\n".format(thousands(GR),
                percentage(GR, G))
        m += "Coverage estimate: {0} x\n".format(coverage)
        self.repetitive = "Repeats: {0} percent".format(GR * 100 / G)

        if ploidy == 2:
            d_SNP = int(_d_SNP)
            self.snprate = "SNP rate ~= 1/{0}".format(d_SNP)
        else:
            self.snprate = "SNP rate not computed (Ploidy = {0})".format(ploidy)
        m += self.snprate + '\n'

        self.genomesize = int(round(self.totalKmers * 1. / self.max2))

        print >> sys.stderr, m


def main():

    actions = (
        ('jellyfish', 'dump histogram using `jellyfish`'),
        ('meryl', 'dump histogram using `meryl`'),
        ('histogram', 'plot the histogram based on meryl K-mer distribution'),
        ('multihistogram', 'plot histogram across a set of K-mer sizes'),
        # These forms a pipeline to count K-mers for given FASTA seq
        ('dump', 'convert FASTA sequences to list of K-mers'),
        ('bin', 'serialize counts to bitarrays'),
        ('bincount', 'count K-mers in the bin'),
        ('count', 'run dump - jellyfish - bin - bincount in serial'),
        ('logodds', 'compute log likelihood between two db'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def logodds(args):
    """
    %prog logodds cnt1 cnt2

    Compute log likelihood between two db.
    """
    from math import log
    from jcvi.formats.base import DictFile

    p = OptionParser(logodds.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    cnt1, cnt2 = args
    d = DictFile(cnt2)
    fp = open(cnt1)
    for row in fp:
        scf, c1 = row.split()
        c2 = d[scf]
        c1, c2 = float(c1), float(c2)
        c1 += 1
        c2 += 1
        score = int(100 * (log(c1) - log(c2)))
        print "{0}\t{1}".format(scf, score)


def get_K(jfdb):
    """
    Infer K from jellyfish db.
    """
    j = jfdb.rsplit('_', 1)[0].rsplit('-', 1)[-1]
    assert j[0] == 'K'
    return int(j[1:])


def count(args):
    """
    %prog count fastafile jf.db

    Run dump - jellyfish - bin - bincount in serial.
    """
    from bitarray import bitarray

    p = OptionParser(count.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, jfdb = args
    K = get_K(jfdb)
    cmd = "jellyfish query {0} -C | cut -d' ' -f 2".format(jfdb)
    t = must_open("tmp", "w")
    proc = Popen(cmd, stdin=PIPE, stdout=t)
    t.flush()

    f = Fasta(fastafile, lazy=True)
    for name, rec in f.iteritems_ordered():
        kmers = list(make_kmers(rec.seq, K))
        print >> proc.stdin, "\n".join(kmers)
    proc.stdin.close()
    logging.debug(cmd)
    proc.wait()

    a = bitarray()
    binfile = ".".join((fastafile, jfdb, "bin"))
    fw = open(binfile, "w")
    t.seek(0)
    for row in t:
        c = row.strip()
        a.append(int(c))
    a.tofile(fw)
    logging.debug("Serialize {0} bits to `{1}`.".format(len(a), binfile))
    fw.close()
    sh("rm {0}".format(t.name))

    logging.debug("Shared K-mers (K={0}) between `{1}` and `{2}` written to `{3}`.".\
                    format(K, fastafile, jfdb, binfile))
    cntfile = ".".join((fastafile, jfdb, "cnt"))
    bincount([fastafile, binfile, "-o", cntfile, "-K {0}".format(K)])
    logging.debug("Shared K-mer counts written to `{0}`.".format(cntfile))


def bincount(args):
    """
    %prog bincount fastafile binfile

    Count K-mers in the bin.
    """
    from bitarray import bitarray
    from jcvi.formats.sizes import Sizes

    p = OptionParser(bincount.__doc__)
    p.add_option("-K", default=23, type="int",
                 help="K-mer size [default: %default]")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, binfile = args
    K = opts.K

    fp = open(binfile)
    a = bitarray()
    a.fromfile(fp)
    f = Sizes(fastafile)
    tsize = 0
    fw = must_open(opts.outfile, "w")
    for name, seqlen in f.iter_sizes():
        ksize = seqlen - K + 1
        b = a[tsize: tsize + ksize]
        bcount = b.count()
        print >> fw, "\t".join(str(x) for x in (name, bcount))
        tsize += ksize


def bin(args):
    """
    %prog bin filename filename.bin

    Serialize counts to bitarrays.
    """
    from bitarray import bitarray
    p = OptionParser(bin.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    inp, outp = args
    fp = must_open(inp)
    fw = must_open(outp, "w")
    a = bitarray()
    for row in fp:
        c = row.split()[-1]
        a.append(int(c))
    a.tofile(fw)
    fw.close()


def make_kmers(seq, K):
    seq = str(seq).upper().replace("N", "A")
    seqlen = len(seq)
    for i in xrange(seqlen - K + 1):
        yield seq[i: i + K]


def dump(args):
    """
    %prog dump fastafile

    Convert FASTA sequences to list of K-mers.
    """
    p = OptionParser(dump.__doc__)
    p.add_option("-K", default=23, type="int",
                 help="K-mer size [default: %default]")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    K = opts.K
    fw = must_open(opts.outfile, "w")
    f = Fasta(fastafile, lazy=True)
    for name, rec in f.iteritems_ordered():
        kmers = list(make_kmers(rec.seq, K))
        print >> fw, "\n".join(kmers)
    fw.close()


def jellyfish(args):
    """
    %prog jellyfish [*.fastq|*.fasta]

    Run jellyfish to dump histogram to be used in kmer.histogram().
    """
    from jcvi.apps.base import getfilesize
    from jcvi.utils.cbook import human_size
    p = OptionParser(jellyfish.__doc__)
    p.add_option("-K", default=23, type="int",
                 help="K-mer size [default: %default]")
    p.add_option("--coverage", default=40, type="int",
                 help="Expected sequence coverage [default: %default]")
    p.add_option("--prefix", default="jf",
                 help="Database prefix [default: %default]")
    p.add_option("--nohist", default=False, action="store_true",
                 help="Do not print histogram [default: %default]")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fastqfiles = args
    K = opts.K
    coverage = opts.coverage

    totalfilesize = sum(getfilesize(x) for x in fastqfiles)
    fq = fastqfiles[0]
    pf = opts.prefix
    gzip = fq.endswith(".gz")

    hashsize = totalfilesize / coverage
    logging.debug("Total file size: {0}, hashsize (-s): {1}".\
                    format(human_size(totalfilesize,
                           a_kilobyte_is_1024_bytes=True), hashsize))

    jfpf = "{0}-K{1}".format(pf, K)
    jfdb = jfpf
    fastqfiles = " ".join(fastqfiles)

    cmd = "jellyfish count -t {0} -C -o {1}".format(opts.cpus, jfpf)
    cmd += " -s {0} -m {1}".format(hashsize, K)
    if gzip:
        cmd = "gzip -dc {0} | ".format(fastqfiles) + cmd + " /dev/fd/0"
    else:
        cmd += " " + fastqfiles

    if need_update(fastqfiles, jfdb):
        sh(cmd)

    if opts.nohist:
        return

    jfhisto = jfpf + ".histogram"
    cmd = "jellyfish histo -t 64 {0} -o {1}".format(jfdb, jfhisto)

    if need_update(jfdb, jfhisto):
        sh(cmd)


def meryl(args):
    """
    %prog meryl merylfile

    Run meryl to dump histogram to be used in kmer.histogram(). The merylfile
    are the files ending in .mcidx or .mcdat.
    """
    p = OptionParser(meryl.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    merylfile, = args
    pf, sf = op.splitext(merylfile)
    outfile = pf + ".histogram"
    cmd = "meryl -Dh -s {0}".format(pf)
    sh(cmd, outfile=outfile)

    return outfile


def multihistogram(args):
    """
    %prog multihistogram *.histogram species

    Plot the histogram based on a set of K-mer hisotograms. The method is based
    on Star et al.'s method (Atlantic Cod genome paper).
    """
    p = OptionParser(multihistogram.__doc__)
    p.add_option("--kmin", default=15, type="int",
            help="Minimum K-mer size, inclusive")
    p.add_option("--kmax", default=30, type="int",
            help="Maximum K-mer size, inclusive")
    p.add_option("--vmin", default=2, type="int",
            help="Minimum value, inclusive")
    p.add_option("--vmax", default=100, type="int",
            help="Maximum value, inclusive")
    opts, args, iopts = p.set_image_options(args, figsize="10x5", dpi=300)

    histfiles = args[:-1]
    species = args[-1]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    A = fig.add_axes([.08, .12, .38, .76])
    B = fig.add_axes([.58, .12, .38, .76])

    lines = []
    legends = []
    genomesizes = []
    for histfile in histfiles:
        ks = KmerSpectrum(histfile)
        x, y = ks.get_xy(opts.vmin, opts.vmax)
        K = get_number(op.basename(histfile).split(".")[0].split("-")[-1])
        if not opts.kmin <= K <= opts.kmax:
            continue

        line, = A.plot(x, y, '-', lw=1)
        lines.append(line)
        legends.append("K = {0}".format(K))
        ks.analyze(K=K)
        genomesizes.append((K, ks.genomesize / 1e6))

    leg = A.legend(lines, legends, shadow=True, fancybox=True)
    leg.get_frame().set_alpha(.5)

    title = "{0} genome K-mer histogram".format(species)
    A.set_title(markup(title))
    xlabel, ylabel = "Coverage (X)", "Counts"
    A.set_xlabel(xlabel)
    A.set_ylabel(ylabel)
    set_human_axis(A)

    title = "{0} genome size estimate".format(species)
    B.set_title(markup(title))
    x, y = zip(*genomesizes)
    B.plot(x, y, "ko", mfc='w')
    t = np.linspace(opts.kmin - .5, opts.kmax + .5, 100)
    p = np.poly1d(np.polyfit(x, y, 2))
    B.plot(t, p(t), "r:")

    xlabel, ylabel = "K-mer size", "Estimated genome size (Mb)"
    B.set_xlabel(xlabel)
    B.set_ylabel(ylabel)
    set_ticklabels_helvetica(B)

    labels = ((.04, .96, 'A'), (.54, .96, 'B'))
    panel_labels(root, labels)

    normalize_axes(root)
    imagename = species + ".multiK.pdf"
    savefig(imagename, dpi=iopts.dpi, iopts=iopts)


def histogram(args):
    """
    %prog histogram meryl.histogram species K

    Plot the histogram based on meryl K-mer distribution, species and N are
    only used to annotate the graphic. Find out totalKmers when running
    kmer.meryl().
    """
    p = OptionParser(histogram.__doc__)
    p.add_option("--vmin", dest="vmin", default=1, type="int",
            help="minimum value, inclusive [default: %default]")
    p.add_option("--vmax", dest="vmax", default=100, type="int",
            help="maximum value, inclusive [default: %default]")
    p.add_option("--pdf", default=False, action="store_true",
            help="Print PDF instead of ASCII plot [default: %default]")
    p.add_option("--coverage", default=0, type="int",
            help="Kmer coverage [default: auto]")
    p.add_option("--nopeaks", default=False, action="store_true",
            help="Do not annotate K-mer peaks")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    histfile, species, N = args
    ascii = not opts.pdf
    peaks = not opts.nopeaks
    N = int(N)

    ks = KmerSpectrum(histfile)
    ks.analyze(K=N)

    Total_Kmers = int(ks.totalKmers)
    coverage = opts.coverage
    Kmer_coverage = ks.max2 if not coverage else coverage
    Genome_size = int(round(Total_Kmers * 1. / Kmer_coverage))

    Total_Kmers_msg = "Total {0}-mers: {1}".format(N, thousands(Total_Kmers))
    Kmer_coverage_msg = "{0}-mer coverage: {1}".format(N, Kmer_coverage)
    Genome_size_msg = "Estimated genome size: {0:.1f}Mb".\
                        format(Genome_size / 1e6)
    Repetitive_msg = ks.repetitive
    SNPrate_msg = ks.snprate

    for msg in (Total_Kmers_msg, Kmer_coverage_msg, Genome_size_msg):
        print >> sys.stderr, msg

    x, y = ks.get_xy(opts.vmin, opts.vmax)
    title = "{0} genome {1}-mer histogram".format(species, N)

    if ascii:
        asciiplot(x, y, title=title)
        return Genome_size

    plt.figure(1, (6, 6))
    plt.plot(x, y, 'g-', lw=2, alpha=.5)
    ax = plt.gca()

    if peaks:
        t = (ks.min1, ks.max1, ks.min2, ks.max2, ks.min3)
        tcounts = [(x, y) for x, y in ks.counts if x in t]
        x, y = zip(*tcounts)
        tcounts = dict(tcounts)
        plt.plot(x, y, 'ko', lw=2, mec='k', mfc='w')
        ax.text(ks.max1, tcounts[ks.max1], "SNP peak", va="top")
        ax.text(ks.max2, tcounts[ks.max2], "Main peak")

    tc = "gray"
    axt = ax.transAxes
    ax.text(.95, .95, Total_Kmers_msg, color=tc, transform=axt, ha="right")
    ax.text(.95, .9, Kmer_coverage_msg, color=tc, transform=axt, ha="right")
    ax.text(.95, .85, Genome_size_msg, color=tc, transform=axt, ha="right")
    ax.text(.95, .8, Repetitive_msg, color=tc, transform=axt, ha="right")
    ax.text(.95, .75, SNPrate_msg, color=tc, transform=axt, ha="right")

    ymin, ymax = ax.get_ylim()
    ymax = ymax * 7 / 6

    ax.set_title(markup(title), color='r')
    ax.set_ylim((ymin, ymax))
    xlabel, ylabel = "Coverage (X)", "Counts"
    ax.set_xlabel(xlabel, color='r')
    ax.set_ylabel(ylabel, color='r')
    set_human_axis(ax)

    imagename = histfile.split(".")[0] + ".pdf"
    savefig(imagename, dpi=100)

    return Genome_size


if __name__ == '__main__':
    main()
