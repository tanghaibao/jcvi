#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Deals with K-mers and K-mer distribution from reads or genome
"""

import os.path as op
import sys
import logging

from optparse import OptionParser

from jcvi.utils.iter import pairwise
from jcvi.graphics.base import plt, asciiplot, _, set_human_axis, savefig
from jcvi.apps.base import ActionDispatcher, sh, debug, need_update
debug()


class KmerSpectrum (object):
    
    def __init__(self, data):
        self.data = data

    def analyze(self, nkf=1000, ploidy=2, K=23):
        """
        Analyze Kmer spectrum, calculations derived from
        allpathslg/src/kmers/KmerSpectra.cc
        """
        from math import sqrt
        from jcvi.utils.cbook import thousands, percentage

        data = self.data
        kf_ceil = nkf - 1
        a = [0] * nkf
        for kf, c in data:
            if kf > kf_ceil:
                break
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

        if _kf_hi > kf_ceil:
            _kf_hi = kf_ceil

        _nk_total = cnk[len(cnk) - 1]
        _nk_bad_low_kf = cnk[_kf_min1]
        _nk_good_snp = cnk[_kf_min2] - cnk[_kf_min1]
        _nk_good_uniq = cnk[_kf_min3] - cnk[_kf_min2]
        _nk_good_rep = cnk[_kf_hi] - cnk[_kf_min3]
        _nk_bad_high_kf = _nk_total - cnk[_kf_hi]

        _ndk_total = cndk[len(cndk) - 1]
        _ndk_bad_low_kf = cndk[_kf_min1]
        _ndk_good_snp = cndk[_kf_min2] - cndk[_kf_min1]
        _ndk_good_uniq = cndk[_kf_min3] - cndk[_kf_min2]
        _ndk_good_rep = cndk[_kf_hi] - cndk[_kf_min3]
        _ndk_bad_high_kf = _ndk_total - cndk[_kf_hi]

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

        m = "Kmer Spectrum Analysis\n"
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

        print >> sys.stderr, m


def main():

    actions = (
        ('jellyfish', 'dump histogram using `jellyfish`'),
        ('meryl', 'dump histogram using `meryl`'),
        ('histogram', 'plot the histogram based on meryl K-mer distribution'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def jellyfish(args):
    """
    %prog jellyfish *.fastq

    Run jellyfish to dump histogram to be used in kmer.histogram().
    """
    from jcvi.apps.base import getfilesize
    from jcvi.utils.cbook import human_size
    from jcvi.formats.fastq import guessoffset

    p = OptionParser(jellyfish.__doc__)
    p.add_option("-K", default=23, type="int",
                 help="K-mer size [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fastqfiles = args
    K = opts.K

    totalfilesize = sum(getfilesize(x) for x in fastqfiles)
    hashsize = totalfilesize / 46
    #hashsize = max(hashsize, 4000000000)  # based on msr-ca

    logging.debug("Total file size: {0}, hashsize (-s): {1}".\
                    format(human_size(totalfilesize,
                           a_kilobyte_is_1024_bytes=True), hashsize))

    offset = guessoffset([fastqfiles[0]])
    assert all(guessoffset([x]) == offset for x in fastqfiles[1:])

    jfpf = "jf-{0}".format(K)
    jfdb = jfpf + "_0"

    cmd = "jellyfish count -t 64 -p 126 -C -o {0}".format(jfpf)
    cmd += " -s {0} -m {1} --min-quality 5".format(hashsize, K)
    cmd += " --quality-start {0}".format(offset)
    cmd += " " + " ".join(fastqfiles)

    if need_update(fastqfiles, jfdb):
        sh(cmd)

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
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    histfile, species, N = args
    N = int(N)
    KMERYL, KSOAP, KALLPATHS = range(3)
    kformats = ("Meryl", "Soap", "AllPaths")
    kformat = KMERYL

    ascii = not opts.pdf
    fp = open(histfile)
    hist = {}
    totalKmers = 0

    # Guess the format of the Kmer histogram
    soap = False
    for row in fp:
        if row.startswith("# 1:"):
            kformat = KALLPATHS
            break
        if len(row.split()) == 1:
            kformat = KSOAP
            break
    fp.seek(0)

    logging.debug("Guessed format: {0}".format(kformats[kformat]))

    data = []
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
        totalKmers += Kcounts
        hist[K] = Kcounts
        data.append((K, counts))

    ks = KmerSpectrum(data)
    ks.analyze(K=N)

    history = ["drop"]
    for a, b in pairwise(sorted(hist.items())):
        Ka, ca = a
        Kb, cb = b
        if ca <= cb:
            status = "rise"
        else:
            status = "drop"
        if history[-1] != status:
            history.append(status)
        if history == ["drop", "rise", "drop"]:
            break

    Total_Kmers = int(totalKmers)
    coverage = opts.coverage
    Kmer_coverage = Ka if not coverage else coverage
    Genome_size = Total_Kmers * 1. / Kmer_coverage / 1e6

    Total_Kmers_msg = "Total {0}-mers: {1}".format(N, Total_Kmers)
    Kmer_coverage_msg = "{0}-mer coverage: {1}".format(N, Kmer_coverage)
    Genome_size_msg = "Estimated genome size: {0:.1f}Mb".format(Genome_size)
    Repetitive_msg = ks.repetitive
    SNPrate_msg = ks.snprate

    for msg in (Total_Kmers_msg, Kmer_coverage_msg, Genome_size_msg):
        print >> sys.stderr, msg

    counts = sorted((a, b) for a, b in hist.items() \
                    if opts.vmin <= a <= opts.vmax)
    x, y = zip(*counts)
    title = "{0} genome {1}-mer histogram".format(species, N)

    if ascii:
        return asciiplot(x, y, title=title)

    fig = plt.figure(1, (6, 6))
    plt.plot(x, y, 'g-', lw=2, alpha=.5)

    t = (ks.min1, ks.max1, ks.min2, ks.max2, ks.min3)
    tcounts = [(x, y) for x, y in counts if x in t]
    x, y = zip(*tcounts)
    plt.plot(x, y, 'ko', lw=2, mec='k', mfc='w')
    tcounts = dict(tcounts)

    ax = plt.gca()

    tc = "gray"
    axt = ax.transAxes
    ax.text(ks.max1, tcounts[ks.max1], "SNP peak", va="top")
    ax.text(ks.max2, tcounts[ks.max2], "Main peak")
    ax.text(.05, .95, _(Total_Kmers_msg), color=tc, transform=axt)
    ax.text(.05, .9, _(Kmer_coverage_msg), color=tc, transform=axt)
    ax.text(.05, .85, _(Genome_size_msg), color=tc, transform=axt)
    ax.text(.05, .8, Repetitive_msg, color=tc, transform=axt)
    ax.text(.05, .75, SNPrate_msg, color=tc, transform=axt)

    ymin, ymax = ax.get_ylim()
    ymax = ymax * 7 / 6

    ax.set_title(_(title), color='r')
    ax.set_ylim((ymin, ymax))
    xlabel, ylabel = "Coverage (X)", "Counts"
    ax.set_xlabel(_(xlabel), color='r')
    ax.set_ylabel(_(ylabel), color='r')
    set_human_axis(ax)

    imagename = histfile.split(".")[0] + ".pdf"
    savefig(imagename, dpi=100)


if __name__ == '__main__':
    main()
