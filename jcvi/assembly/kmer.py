#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Deals with K-mers and K-mer distribution from reads or genome
"""
import os.path as op
import sys
import math

from collections import defaultdict
from typing import List

import numpy as np
from more_itertools import chunked

from ..apps.grid import MakeManager
from ..apps.base import (
    ActionDispatcher,
    OptionParser,
    PIPE,
    Popen,
    logger,
    need_update,
    sh,
)
from ..formats.fasta import Fasta
from ..formats.base import BaseFile, must_open, get_number
from ..graphics.base import (
    adjust_spines,
    asciiplot,
    markup,
    normalize_axes,
    panel_labels,
    plt,
    savefig,
    set_human_axis,
    set_ticklabels_helvetica,
    write_messages,
)
from ..utils.cbook import thousands, percentage

from .automaton import iter_project


KMERYL, KSOAP, KALLPATHS = range(3)


class KmerSpectrum(BaseFile):
    def __init__(self, histfile):
        super(KmerSpectrum, self).__init__(histfile)
        self.load_data(histfile)

    def load_data(self, histfile):
        self.data = []
        self.totalKmers = 0
        self.hist = {}
        kformat = self.guess_format(histfile)
        kformats = ("Meryl", "Soap", "AllPaths")
        logger.debug("Guessed format: %s", kformats[kformat])

        fp = open(histfile)
        for rowno, row in enumerate(fp):
            if row[0] == "#":
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
        self.counts = sorted((a, b) for a, b in self.hist.items() if vmin <= a <= vmax)
        return zip(*self.counts)

    def analyze(self, K=23, maxiter=100, method="nbinom"):
        """Analyze K-mer histogram.

        Args:
            K (int, optional): K-mer size. Defaults to 23.
            maxiter (int): Iterations to run. Defaults to 100.
            method (str, optional): Method to use, either 'nbinom' or
            'allpaths'. Defaults to "nbinom".

        Returns:
            A dictionary containing info for annotating the plot. analyze() also
            sets the following properties:
            - lambda_: Main peak
            - repetitive: Repeats message
            - snprate: SNP rate message
        """
        if method == "nbinom":
            return self.analyze_nbinom(K=K, maxiter=maxiter)
        return self.analyze_allpaths(K=K)

    def analyze_nbinom(self, K=23, maxiter=100):
        """Analyze the K-mer histogram using negative binomial distribution.

        Args:
            K (int, optional): K-mer size used when generating the histogram. Defaults to 23.
        """
        from scipy.stats import nbinom
        from scipy.optimize import minimize_scalar
        from functools import lru_cache

        method, xopt = "bounded", "xatol"
        MAX_1CN_SIZE = 1e10
        MAX_OPTIMIZED_SIZE = 9.9e9

        # Generate bins for the decomposed negative binomial distributions
        bins = [
            (i, i) for i in range(1, 9)
        ]  # The first 8 CN are critical often determines ploidy
        for i in (8, 16, 32, 64, 128, 256, 512):  # 14 geometricly sized bins
            a, b = i + 1, int(round(i * 2**0.5))
            bins.append((a, b))
            a, b = b + 1, i * 2
            bins.append((a, b))

        # Convert histogram to np array so we can index by CN
        kf_ceil = max([cov for cov, _ in self.data])
        N = kf_ceil + 1
        hist = np.zeros(N, dtype=int)
        for cov, count in self.data:
            hist[cov] = count

        # min1: find first minimum
        _kf_min1 = 5
        while (
            _kf_min1 - 1 >= 2
            and hist[_kf_min1 - 1] * (_kf_min1 - 1) < hist[_kf_min1] * _kf_min1
        ):
            _kf_min1 -= 1
        while (
            _kf_min1 <= kf_ceil
            and hist[_kf_min1 + 1] * (_kf_min1 + 1) < hist[_kf_min1] * _kf_min1
        ):
            _kf_min1 += 1

        # max2: find absolute maximum mx2 above first minimum min1
        _kf_max2 = _kf_min1
        for kf in range(_kf_min1 + 1, int(0.8 * kf_ceil)):
            if hist[kf] * kf > hist[_kf_max2] * _kf_max2:
                _kf_max2 = kf

        # Discard the last entry as that is usually an inflated number
        hist = hist[:-1]
        kf_range = np.arange(_kf_min1, len(hist), dtype=int)
        P = hist[kf_range] * kf_range  # Target distribution
        print("==> Start nbinom method on range ({}, {})".format(_kf_min1, len(hist)))

        # Below is the optimization schemes, we optimize one variable at a time
        @lru_cache(maxsize=None)
        def nbinom_pmf_range(lambda_: int, rho: int, bin_id: int):
            stacked = np.zeros(len(kf_range), dtype=np.float64)
            lambda_ /= 100  # 2-digit precision
            rho /= 100  # 2-digit precision
            n = lambda_ / (rho - 1)
            p = 1 / rho
            start, end = bins[bin_id]
            for i in range(start, end + 1):
                stacked += nbinom.pmf(kf_range, n * i, p)
            return stacked

        def generative_model(G, lambda_, rho):
            stacked = np.zeros(len(kf_range), dtype=np.float64)
            lambda_ = int(round(lambda_ * 100))
            rho = int(round(rho * 100))
            for bin_id, g in enumerate(G):
                stacked += g * nbinom_pmf_range(lambda_, rho, bin_id)
            stacked *= kf_range
            return stacked

        def func(lambda_, rho, G):
            stacked = generative_model(G, lambda_, rho)
            return np.sum((P - stacked) ** 2)  # L2 norm

        def optimize_func(lambda_, rho, G):
            # Iterate over all G
            for i, g in enumerate(G):
                G_i = optimize_func_Gi(lambda_, rho, G, i)
                if (
                    not 1 < G_i < MAX_OPTIMIZED_SIZE
                ):  # Optimizer did not optimize this G_i
                    break
            # Also remove the last bin since it is subject to marginal effect
            G[i - 1] = 0
            lambda_ = optimize_func_lambda_(lambda_, rho, G)
            rho = optimize_func_rho(lambda_, rho, G)
            score = func(lambda_, rho, G)
            return lambda_, rho, G, score

        def optimize_func_lambda_(lambda_, rho, G):
            def f(arg):
                return func(arg, rho, G)

            res = minimize_scalar(
                f, bounds=(_kf_min1, 100), method=method, options={xopt: 0.01}
            )
            return res.x

        def optimize_func_rho(lambda_, rho, G):
            def f(arg):
                return func(lambda_, arg, G)

            res = minimize_scalar(
                f, bounds=(1.001, 5), method=method, options={xopt: 0.01}
            )
            return res.x

        def optimize_func_Gi(lambda_, rho, G, i):
            # Iterate a single G_i
            def f(arg):
                G[i] = arg
                return func(lambda_, rho, G)

            res = minimize_scalar(
                f, bounds=(0, MAX_1CN_SIZE), method=method, options={xopt: 100}
            )
            return res.x

        def run_optimization(termination=0.999, maxiter=100):
            ll, rr, GG = l0, r0, G0
            prev_score = np.inf
            for i in range(maxiter):
                print("Iteration", i + 1, file=sys.stderr)
                ll, rr, GG, score = optimize_func(ll, rr, GG)
                if score / prev_score > termination:
                    break
                prev_score = score
                if i % 10 == 0:
                    print(ll, rr, GG, score, file=sys.stderr)
            print("Success!", file=sys.stderr)
            # Remove bogus values that are close to the bounds
            final_GG = [g for g in GG if 1 < g < MAX_OPTIMIZED_SIZE]
            return ll, rr, final_GG

        # Optimization - very slow
        G0 = np.zeros(len(bins))
        l0 = _kf_max2
        r0 = 1.5
        print(l0, r0, G0, file=sys.stderr)
        ll, rr, GG = run_optimization(maxiter=maxiter)
        print(ll, rr, GG, file=sys.stderr)

        # Ready for genome summary
        m = "\n==> Kmer (K={0}) Spectrum Analysis\n".format(K)

        genome_size = int(round(self.totalKmers / ll))
        inferred_genome_size = 0
        for i, g in enumerate(GG):
            start, end = bins[i]
            mid = (start + end) / 2
            inferred_genome_size += g * mid * (end - start + 1)
        inferred_genome_size = int(round(inferred_genome_size))
        genome_size = max(genome_size, inferred_genome_size)
        m += "Genome size estimate = {0}\n".format(thousands(genome_size))
        copy_series = []
        copy_messages = []
        for i, g in enumerate(GG):
            start, end = bins[i]
            mid = (start + end) / 2
            copy_num = start if start == end else "{}-{}".format(start, end)
            g_copies = int(round(g * mid * (end - start + 1)))
            copy_series.append((mid, copy_num, g_copies, g))
            copy_message = f"CN {copy_num}: {g_copies / 1e6:.1f} Mb ({ g_copies * 100 / genome_size:.1f} percent)"
            copy_messages.append(copy_message)
            m += copy_message + "\n"

        if genome_size > inferred_genome_size:
            g_copies = genome_size - inferred_genome_size
            copy_num = "{}+".format(end + 1)
            copy_series.append((end + 1, copy_num, g_copies, g_copies / (end + 1)))
            m += f"CN {copy_num}: {g_copies / 1e6:.1f} Mb ({ g_copies * 100 / genome_size:.1f} percent)\n"

        # Determine ploidy
        def determine_ploidy(copy_series, threshold=0.15):
            counts_so_far = 1
            ploidy_so_far = 0
            for mid, _, g_copies, _ in copy_series:
                if g_copies / counts_so_far < threshold:
                    break
                counts_so_far += g_copies
                ploidy_so_far = mid
            return int(ploidy_so_far)

        ploidy = determine_ploidy(copy_series)
        self.ploidy = ploidy
        self.ploidy_message = "Ploidy: {}".format(ploidy)
        m += self.ploidy_message + "\n"
        self.copy_messages = copy_messages[:ploidy]

        # Repeat content
        def calc_repeats(copy_series, ploidy, genome_size):
            unique = 0
            for mid, _, g_copies, _ in copy_series:
                if mid <= ploidy:
                    unique += g_copies
                else:
                    break
            return 1 - unique / genome_size

        repeats = calc_repeats(copy_series, ploidy, genome_size)
        self.repetitive = "Repeats: {:.1f} percent".format(repeats * 100)
        m += self.repetitive + "\n"

        # SNP rate
        def calc_snp_rate(copy_series, ploidy, genome_size, K):
            # We can calculate the SNP rate s, assuming K-mer of length K:
            # s = 1-(1-L/G)^(1/K)
            # L: # of unique K-mers under 'het' peak
            # G: genome size
            # K: K-mer length
            L = 0
            for mid, copy_num, g_copies, g in copy_series:
                if mid < ploidy:
                    L += g
                else:
                    break
            return 1 - (1 - L / genome_size) ** (1 / K)

        snp_rate = calc_snp_rate(copy_series, ploidy, genome_size, K)
        self.snprate = "SNP rate: {:.2f} percent".format(snp_rate * 100)
        m += self.snprate + "\n"
        print(m, file=sys.stderr)

        self.lambda_ = ll
        return {
            "generative_model": generative_model,
            "Gbins": GG,
            "lambda": ll,
            "rho": rr,
            "kf_range": kf_range,
        }

    def analyze_allpaths(self, ploidy=2, K=23, covmax=1000000):
        """
        Analyze Kmer spectrum, calculations derived from
        allpathslg/src/kmers/KmerSpectra.cc
        """
        from math import sqrt

        data = self.data
        kf_ceil = max(K for (K, c) in data)
        if kf_ceil > covmax:
            exceeds = sum(1 for (K, c) in data if K > covmax)
            logger.debug(
                "A total of %d distinct K-mers appear > %d times. Ignored ...",
                exceeds,
                covmax,
            )
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
        for kf in range(1, nkf):
            cndk[kf] = cndk[kf - 1] + 0.5 * (ndk[kf - 1] + ndk[kf])
            cnk[kf] = cnk[kf - 1] + 0.5 * (nk[kf - 1] + nk[kf])

        # Separate kmer spectrum in 5 regions based on the kf
        # 1        ... kf_min1    : bad kmers with low frequency
        # kf_min1  ... kf_min2    : good kmers CN = 1/2 (SNPs)
        # kf_min2  ... kf_min3    : good kmers CN = 1
        # kf_min3  ... kf_hi      : good kmers CN > 1 (repetitive)
        # kf_hi    ... inf        : bad kmers with high frequency

        # min1: find first minimum
        _kf_min1 = 10
        while _kf_min1 - 1 >= 2 and nk[_kf_min1 - 1] < nk[_kf_min1]:
            _kf_min1 -= 1
        while _kf_min1 <= kf_ceil and nk[_kf_min1 + 1] < nk[_kf_min1]:
            _kf_min1 += 1

        # max2: find absolute maximum mx2 above first minimum min1
        _kf_max2 = _kf_min1
        for kf in range(_kf_min1 + 1, int(0.8 * kf_ceil)):
            if nk[kf] > nk[_kf_max2]:
                _kf_max2 = kf

        # max2: resetting max2 for cases of very high polymorphism
        if ploidy == 2:
            ndk_half = ndk[_kf_max2 // 2]
            ndk_double = ndk[_kf_max2 * 2]
            if ndk_double > ndk_half:
                _kf_max2 *= 2

        # max1: SNPs local maximum max1 as half global maximum max2
        _kf_max1 = _kf_max2 // 2

        # min2: SNPs local minimum min2 between max1 and max2
        _kf_min2 = (
            _kf_max1
            * (2 * ndk[_kf_max1] + ndk[_kf_max2])
            // (ndk[_kf_max1] + ndk[_kf_max2])
        )

        # min1: refine between min1 and max2/2
        for kf in range(_kf_min1 + 1, _kf_max1):
            if nk[kf] < nk[_kf_min1]:
                _kf_min1 = kf

        # min3: not a minimum, really. upper edge of main peak
        _kf_min3 = _kf_max2 * 3 // 2

        print("kfs:", _kf_min1, _kf_max1, _kf_min2, _kf_max2, _kf_min3, file=sys.stderr)
        self.min1 = _kf_min1
        self.max1 = _kf_max1
        self.min2 = _kf_min2
        self.max2 = _kf_max2
        self.min3 = _kf_min3
        self.lambda_ = self.max2  # Main peak

        # Define maximum kf above which we neglect data
        _kf_hi = (
            _kf_max2 * sqrt(4 * ndk[2 * _kf_max2] * _kf_max2)
            if 2 * _kf_max2 < len(ndk)
            else _kf_max2 * sqrt(4 * ndk[len(ndk) - 1] * _kf_max2)
        )
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
        _kf_ave_uniq = _nk_good_uniq * 1.0 / _ndk_good_uniq
        _genome_size = (_nk_total - _nk_bad_low_kf - _nk_bad_high_kf) / _kf_ave_uniq
        _genome_size_unique = _ndk_good_uniq + _ndk_good_snp / 2
        _genome_size_repetitive = _genome_size - _genome_size_unique
        _coverage = _nk_total / _genome_size if _genome_size else 0

        # SNP rate estimation, assumes uniform distribution of SNPs over the
        # genome and accounts for the reduction in SNP kmer counts when
        # polymorphism is very high
        if ploidy == 2:
            _d_SNP = (
                1.0 / (1.0 - (1.0 - 0.5 * _ndk_good_snp / _genome_size) ** (1.0 / K))
                if _ndk_good_snp > 0
                else 1000000
            )

        G = int(_genome_size)
        G1 = int(_genome_size_unique)
        GR = int(_genome_size_repetitive)
        coverage = int(_coverage)

        m = "Kmer (K={0}) Spectrum Analysis\n".format(K)
        m += "Genome size estimate = {0}\n".format(thousands(G))
        m += "Genome size estimate CN = 1 = {0} ({1})\n".format(
            thousands(G1), percentage(G1, G)
        )
        m += "Genome size estimate CN > 1 = {0} ({1})\n".format(
            thousands(GR), percentage(GR, G)
        )
        m += "Coverage estimate: {0} x\n".format(coverage)
        self.repetitive = "Repeats: {0} percent".format(GR * 100 // G)

        if ploidy == 2:
            d_SNP = int(_d_SNP)
            self.snprate = "SNP rate ~= 1/{0}".format(d_SNP)
        else:
            self.snprate = "SNP rate not computed (Ploidy = {0})".format(ploidy)
        m += self.snprate + "\n"

        self.genomesize = int(round(self.totalKmers * 1.0 / self.max2))

        print(m, file=sys.stderr)
        return {}


class KMCComplex(object):
    def __init__(self, indices):
        self.indices = indices

    def write(
        self,
        outfile: str,
        action: str = "union",
        ci_in: int = 0,
        ci_out: int = 0,
        batch: int = 0,
    ):
        assert action in ("union", "intersect")
        op = " + sum " if action == "union" else " * "
        mm = MakeManager()
        if batch > 1:
            filename = outfile + ".{}.def"
            # Divide indices into batches
            batches = []
            batchsize = (len(self.indices) + batch - 1) // batch
            logger.debug("Use batchsize of %d", batchsize)
            for i, indices in enumerate(chunked(self.indices, batchsize)):
                filename_i = filename.format(i + 1)
                outfile_i = outfile + ".{}".format(i + 1)
                self.write_definitions(
                    filename_i, indices, outfile_i, op, ci_in=ci_in, ci_out=0
                )
                cmd = "kmc_tools complex {}".format(filename_i)
                outfile_suf = outfile_i + ".kmc_suf"
                mm.add(indices, outfile_suf, cmd)
                batches.append(outfile_suf)
        else:
            batches = self.indices

        # Merge batches into one
        filename = outfile + ".def"
        self.write_definitions(
            filename, batches, outfile, op, ci_in=ci_in, ci_out=ci_out
        )
        outfile_suf = outfile + ".kmc_suf"
        mm.add(batches, outfile_suf, "kmc_tools complex {}".format(filename))

        # Write makefile
        mm.write()

    def write_definitions(
        self,
        filename: str,
        indices: List[str],
        outfile: str,
        op: str,
        ci_in: int,
        ci_out: int,
    ):
        fw = must_open(filename, "w")
        print("INPUT:", file=fw)
        ss = []
        pad = len(str(len(indices)))
        for i, e in enumerate(indices):
            s = "s{0:0{1}d}".format(i + 1, pad)
            ss.append(s)
            msg = "{} = {}".format(s, e.rsplit(".", 1)[0])
            if ci_in:
                msg += f" -ci{ci_in}"
            print(msg, file=fw)
        print("OUTPUT:", file=fw)
        print("{} = {}".format(outfile, op.join(ss)), file=fw)
        if ci_out:
            print("OUTPUT_PARAMS:", file=fw)
            print(f"-ci{ci_out}", file=fw)
        fw.close()


def main():

    actions = (
        # K-mer counting
        ("jellyfish", "count kmers using `jellyfish`"),
        ("meryl", "count kmers using `meryl`"),
        ("kmc", "count kmers using `kmc`"),
        ("kmcop", "intersect or union kmc indices"),
        ("entropy", "calculate entropy for kmers from kmc dump"),
        ("bed", "map kmers on FASTA"),
        # K-mer histogram
        ("histogram", "plot the histogram based on meryl K-mer distribution"),
        ("multihistogram", "plot histogram across a set of K-mer sizes"),
        # These forms a pipeline to count K-mers for given FASTA seq
        ("dump", "convert FASTA sequences to list of K-mers"),
        ("bin", "serialize counts to bitarrays"),
        ("bincount", "count K-mers in the bin"),
        ("count", "run dump - jellyfish - bin - bincount in serial"),
        ("logodds", "compute log likelihood between two db"),
        ("model", "model kmer distribution given error rate"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def entropy_score(kmer):
    """
    Schmieder and Edwards. Quality control and preprocessing of metagenomic datasets. (2011) Bioinformatics
    https://academic.oup.com/bioinformatics/article/27/6/863/236283/Quality-control-and-preprocessing-of-metagenomic
    """
    l = len(kmer) - 2
    k = l if l < 64 else 64
    counts = defaultdict(int)
    for i in range(l):
        trinuc = kmer[i : i + 3]
        counts[trinuc] += 1

    logk = math.log(k)
    res = 0
    for k, v in counts.items():
        f = v * 1.0 / l
        res += f * math.log(f) / logk
    return res * -100


def entropy(args):
    """
    %prog entropy kmc_dump.out

    kmc_dump.out contains two columns:
    AAAAAAAAAAAGAAGAAAGAAA  34
    """
    p = OptionParser(entropy.__doc__)
    p.add_argument(
        "--threshold", default=0, type=int, help="Complexity needs to be above"
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (kmc_out,) = args
    fp = open(kmc_out)
    for row in fp:
        kmer, count = row.split()
        score = entropy_score(kmer)
        if score >= opts.threshold:
            print(" ".join((kmer, count, "{:.2f}".format(score))))


def bed(args):
    """
    %prog bed fastafile kmer.dump.txt

    Map kmers on FASTA.
    """
    from jcvi.formats.fasta import rc, parse_fasta

    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, dumpfile = args
    fp = open(dumpfile)
    KMERS = set()
    for row in fp:
        kmer = row.split()[0]
        kmer_rc = rc(kmer)
        KMERS.add(kmer)
        KMERS.add(kmer_rc)

    K = len(kmer)
    logger.debug("Imported %d %d-mers", len(KMERS), K)

    for name, seq in parse_fasta(fastafile):
        name = name.split()[0]
        for i in range(len(seq) - K):
            if i % 5000000 == 0:
                print("{}:{}".format(name, i), file=sys.stderr)
            kmer = seq[i : i + K]
            if kmer in KMERS:
                print("\t".join(str(x) for x in (name, i, i + K, kmer)))


def kmcop(args):
    """
    %prog kmcop *.kmc_suf

    Intersect or union kmc indices.
    """
    p = OptionParser(kmcop.__doc__)
    p.add_argument(
        "--action",
        choices=("union", "intersect", "reduce"),
        default="union",
        help="Action",
    )
    p.add_argument(
        "--ci_in",
        default=0,
        type=int,
        help="Exclude input kmers with less than ci_in counts",
    )
    p.add_argument(
        "--cs",
        default=0,
        type=int,
        help="Maximal value of a counter, only used when action is reduce",
    )
    p.add_argument(
        "--ci_out",
        default=0,
        type=int,
        help="Exclude output kmers with less than ci_out counts",
    )
    p.add_argument(
        "--batch",
        default=1,
        type=int,
        help="Number of batch, useful to reduce memory usage",
    )
    p.add_argument("--exclude", help="Exclude accessions from this list")
    p.add_argument("-o", default="results", help="Output name")
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    indices = args
    if opts.exclude:
        before = set(indices)
        exclude_ids = set(x.strip() for x in open(opts.exclude))
        indices = [x for x in indices if x.rsplit(".", 2)[0] not in exclude_ids]
        after = set(indices)
        if before > after:
            logger.debug(
                "Excluded accessions %d → %d (%s)",
                len(before),
                len(after),
                ",".join(before - after),
            )
    if opts.action == "reduce":
        mm = MakeManager()
        ci = opts.ci_in
        cs = opts.cs
        suf = ""
        if ci:
            suf += f"_ci{ci}"
        if cs:
            suf += f"_cs{cs}"
        for index in indices:
            idx = index.rsplit(".", 1)[0]
            reduced_idx = idx + suf
            cmd = f"kmc_tools transform {idx} reduce {reduced_idx}"
            if ci:
                cmd += f" -ci{ci}"
            if cs:
                cmd += f" -cs{cs}"
            reduced_index = reduced_idx + ".kmc_suf"
            mm.add(index, reduced_index, cmd)
        mm.write()
    else:
        ku = KMCComplex(indices)
        ku.write(
            opts.o,
            action=opts.action,
            ci_in=opts.ci_in,
            ci_out=opts.ci_out,
            batch=opts.batch,
        )


def kmc(args):
    """
    %prog kmc folder

    Run kmc3 on Illumina reads.
    """
    p = OptionParser(kmc.__doc__)
    p.add_argument("-k", default=27, type=int, help="Kmer size")
    p.add_argument(
        "--ci", default=2, type=int, help="Exclude kmers with less than ci counts"
    )
    p.add_argument("--cs", default=0, type=int, help="Maximal value of a counter")
    p.add_argument("--cx", type=int, help="Exclude kmers with more than cx counts")
    p.add_argument(
        "--single",
        default=False,
        action="store_true",
        help="Input is single-end data, only one FASTQ/FASTA",
    )
    p.add_argument(
        "--fasta",
        default=False,
        action="store_true",
        help="Input is FASTA instead of FASTQ",
    )
    p.add_argument(
        "--mem", default=48, type=int, help="Max amount of RAM in GB (`kmc -m`)"
    )
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (folder,) = args
    K = opts.k
    n = 1 if opts.single else 2
    pattern = (
        "*.fa,*.fa.gz,*.fasta,*.fasta.gz"
        if opts.fasta
        else "*.fq,*.fq.gz,*.fastq,*.fastq.gz"
    )

    mm = MakeManager()
    for p, pf in iter_project(folder, pattern=pattern, n=n, commonprefix=False):
        pf = pf.split("_")[0] + ".ms{}".format(K)
        infiles = pf + ".infiles"
        fw = open(infiles, "w")
        print("\n".join(p), file=fw)
        fw.close()

        cmd = "kmc -k{} -m{} -t{}".format(K, opts.mem, opts.cpus)
        cmd += " -ci{}".format(opts.ci)
        if opts.cs:
            cmd += " -cs{}".format(opts.cs)
        if opts.cx:
            cmd += " -cx{}".format(opts.cx)
        if opts.fasta:
            cmd += " -fm"
        cmd += " @{} {} .".format(infiles, pf)
        outfile = pf + ".kmc_suf"
        mm.add(p, outfile, cmd)

    mm.write()


def meryl(args):
    """
    %prog meryl folder

    Run meryl on Illumina reads.
    """
    p = OptionParser(meryl.__doc__)
    p.add_argument("-k", default=19, type=int, help="Kmer size")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (folder,) = args
    K = opts.k
    cpus = opts.cpus
    mm = MakeManager()
    for p, pf in iter_project(folder):
        cmds = []
        mss = []
        for i, ip in enumerate(p):
            ms = "{}{}.ms{}".format(pf, i + 1, K)
            mss.append(ms)
            cmd = "meryl -B -C -m {} -threads {}".format(K, cpus)
            cmd += " -s {} -o {}".format(ip, ms)
            cmds.append(cmd)
        ams, bms = mss
        pms = "{}.ms{}".format(pf, K)
        cmd = "meryl -M add -s {} -s {} -o {}".format(ams, bms, pms)
        cmds.append(cmd)
        cmd = "rm -f {}.mcdat {}.mcidx {}.mcdat {}.mcidx".format(ams, ams, bms, bms)
        cmds.append(cmd)
        mm.add(p, pms + ".mcdat", cmds)

    mm.write()


def model(args):
    """
    %prog model erate

    Model kmer distribution given error rate. See derivation in FIONA paper:
    <http://bioinformatics.oxfordjournals.org/content/30/17/i356.full>
    """
    from scipy.stats import binom, poisson

    p = OptionParser(model.__doc__)
    p.add_argument("-k", default=23, type=int, help="Kmer size")
    p.add_argument("--cov", default=50, type=int, help="Expected coverage")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (erate,) = args
    erate = float(erate)
    cov = opts.cov
    k = opts.k

    xy = []
    # Range include c although it is unclear what it means to have c=0
    for c in range(0, cov * 2 + 1):
        Prob_Yk = 0
        for i in range(k + 1):
            # Probability of having exactly i errors
            pi_i = binom.pmf(i, k, erate)
            # Expected coverage of kmer with exactly i errors
            mu_i = cov * (erate / 3) ** i * (1 - erate) ** (k - i)
            # Probability of seeing coverage of c
            Prob_Yk_i = poisson.pmf(c, mu_i)
            # Sum i over 0, 1, ... up to k errors
            Prob_Yk += pi_i * Prob_Yk_i
        xy.append((c, Prob_Yk))

    x, y = zip(*xy)
    asciiplot(x, y, title="Model")


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
        print("{0}\t{1}".format(scf, score))


def get_K(jfdb):
    """
    Infer K from jellyfish db.
    """
    j = jfdb.rsplit("_", 1)[0].rsplit("-", 1)[-1]
    assert j[0] == "K"
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
        print("\n".join(kmers), file=proc.stdin)
    proc.stdin.close()
    logger.debug(cmd)
    proc.wait()

    a = bitarray()
    binfile = ".".join((fastafile, jfdb, "bin"))
    fw = open(binfile, "w")
    t.seek(0)
    for row in t:
        c = row.strip()
        a.append(int(c))
    a.tofile(fw)
    logger.debug("Serialize %d bits to `%s`.", len(a), binfile)
    fw.close()
    sh("rm {0}".format(t.name))

    logger.debug(
        "Shared K-mers (K=%d) between `%s` and `%s` written to `%s`.",
        K,
        fastafile,
        jfdb,
        binfile,
    )
    cntfile = ".".join((fastafile, jfdb, "cnt"))
    bincount([fastafile, binfile, "-o", cntfile, "-K {0}".format(K)])
    logger.debug("Shared K-mer counts written to `%s`.", cntfile)


def bincount(args):
    """
    %prog bincount fastafile binfile

    Count K-mers in the bin.
    """
    from bitarray import bitarray
    from jcvi.formats.sizes import Sizes

    p = OptionParser(bincount.__doc__)
    p.add_argument("-K", default=23, type=int, help="K-mer size")
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
        b = a[tsize : tsize + ksize]
        bcount = b.count()
        print("\t".join(str(x) for x in (name, bcount)), file=fw)
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
    for i in range(seqlen - K + 1):
        yield seq[i : i + K]


def dump(args):
    """
    %prog dump fastafile

    Convert FASTA sequences to list of K-mers.
    """
    p = OptionParser(dump.__doc__)
    p.add_argument("-K", default=23, type=int, help="K-mer size")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    K = opts.K
    fw = must_open(opts.outfile, "w")
    f = Fasta(fastafile, lazy=True)
    for name, rec in f.iteritems_ordered():
        kmers = list(make_kmers(rec.seq, K))
        print("\n".join(kmers), file=fw)
    fw.close()


def jellyfish(args):
    """
    %prog jellyfish [*.fastq|*.fasta]

    Run jellyfish to dump histogram to be used in kmer.histogram().
    """
    from jcvi.apps.base import getfilesize
    from jcvi.utils.cbook import human_size

    p = OptionParser(jellyfish.__doc__)
    p.add_argument("-K", default=23, type=int, help="K-mer size")
    p.add_argument(
        "--coverage",
        default=40,
        type=int,
        help="Expected sequence coverage",
    )
    p.add_argument("--prefix", default="jf", help="Database prefix")
    p.add_argument(
        "--nohist",
        default=False,
        action="store_true",
        help="Do not print histogram",
    )
    p.set_home("jellyfish")
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
    logger.debug(
        "Total file size: %s, hashsize (-s): %d",
        human_size(totalfilesize, a_kilobyte_is_1024_bytes=True),
        hashsize,
    )

    jfpf = "{0}-K{1}".format(pf, K)
    jfdb = jfpf
    fastqfiles = " ".join(fastqfiles)

    jfcmd = op.join(opts.jellyfish_home, "jellyfish")
    cmd = jfcmd
    cmd += " count -t {0} -C -o {1}".format(opts.cpus, jfpf)
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
    cmd = jfcmd + " histo -t 64 {0} -o {1}".format(jfdb, jfhisto)

    if need_update(jfdb, jfhisto):
        sh(cmd)


def multihistogram(args):
    """
    %prog multihistogram *.histogram species

    Plot the histogram based on a set of K-mer hisotograms. The method is based
    on Star et al.'s method (Atlantic Cod genome paper).
    """
    p = OptionParser(multihistogram.__doc__)
    p.add_argument("--kmin", default=15, type=int, help="Minimum K-mer size, inclusive")
    p.add_argument("--kmax", default=30, type=int, help="Maximum K-mer size, inclusive")
    p.add_argument("--vmin", default=2, type=int, help="Minimum value, inclusive")
    p.add_argument("--vmax", default=100, type=int, help="Maximum value, inclusive")
    opts, args, iopts = p.set_image_options(args, figsize="10x5", dpi=300)

    if len(args) < 1:
        sys.exit(not p.print_help())

    histfiles = args[:-1]
    species = args[-1]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))
    A = fig.add_axes((0.08, 0.12, 0.38, 0.76))
    B = fig.add_axes((0.58, 0.12, 0.38, 0.76))

    lines = []
    legends = []
    genomesizes = []
    for histfile in histfiles:
        ks = KmerSpectrum(histfile)
        x, y = ks.get_xy(opts.vmin, opts.vmax)
        K = get_number(op.basename(histfile).split(".")[0].split("-")[-1])
        if not opts.kmin <= K <= opts.kmax:
            continue

        (line,) = A.plot(x, y, "-", lw=1)
        lines.append(line)
        legends.append("K = {0}".format(K))
        ks.analyze(K=K, method="allpaths")
        genomesizes.append((K, ks.genomesize / 1e6))

    leg = A.legend(lines, legends, shadow=True, fancybox=True)
    leg.get_frame().set_alpha(0.5)

    title = "{0} genome K-mer histogram".format(species)
    A.set_title(markup(title))
    xlabel, ylabel = "Coverage (X)", "Counts"
    A.set_xlabel(xlabel)
    A.set_ylabel(ylabel)
    set_human_axis(A)

    title = "{0} genome size estimate".format(species)
    B.set_title(markup(title))
    x, y = zip(*genomesizes)
    B.plot(x, y, "ko", mfc="w")
    t = np.linspace(opts.kmin - 0.5, opts.kmax + 0.5, 100)
    p = np.poly1d(np.polyfit(x, y, 2))
    B.plot(t, p(t), "r:")

    xlabel, ylabel = "K-mer size", "Estimated genome size (Mb)"
    B.set_xlabel(xlabel)
    B.set_ylabel(ylabel)
    set_ticklabels_helvetica(B)

    labels = ((0.04, 0.96, "A"), (0.54, 0.96, "B"))
    panel_labels(root, labels)

    normalize_axes(root)
    imagename = species + ".multiK.pdf"
    savefig(imagename, dpi=iopts.dpi, iopts=iopts)


def plot_nbinom_fit(ax, ks: KmerSpectrum, ymax: float, method_info: dict):
    """
    Plot the negative binomial fit.
    """
    generative_model = method_info["generative_model"]
    GG = method_info["Gbins"]
    ll = method_info["lambda"]
    rr = method_info["rho"]
    kf_range = method_info["kf_range"]
    stacked = generative_model(GG, ll, rr)
    ax.plot(
        kf_range,
        stacked,
        ":",
        color="#6a3d9a",
        lw=2,
    )
    # Plot multiple CN locations, CN1, CN2, ... up to ploidy
    cn_color = "#a6cee3"
    for i in range(1, ks.ploidy + 1):
        x = i * ks.lambda_
        ax.plot((x, x), (0, ymax), "-.", color=cn_color)
        ax.text(
            x,
            ymax * 0.95,
            f"CN{i}",
            ha="right",
            va="center",
            color=cn_color,
            rotation=90,
        )


def draw_ks_histogram(
    ax,
    histfile: str,
    method: str,
    coverage: int,
    vmin: int,
    vmax: int,
    species: str,
    K: int,
    maxiter: int,
    peaks: bool,
) -> int:
    """
    Draw the K-mer histogram.
    """
    ks = KmerSpectrum(histfile)
    method_info = ks.analyze(K=K, maxiter=maxiter, method=method)

    Total_Kmers = int(ks.totalKmers)
    Kmer_coverage = ks.lambda_ if not coverage else coverage
    Genome_size = int(round(Total_Kmers * 1.0 / Kmer_coverage))

    Total_Kmers_msg = f"Total {K}-mers: {thousands(Total_Kmers)}"
    Kmer_coverage_msg = f"{K}-mer coverage: {Kmer_coverage:.1f}x"
    Genome_size_msg = f"Estimated genome size: {Genome_size / 1e6:.1f} Mb"
    Repetitive_msg = ks.repetitive
    SNPrate_msg = ks.snprate

    messages = [
        Total_Kmers_msg,
        Kmer_coverage_msg,
        Genome_size_msg,
        Repetitive_msg,
        SNPrate_msg,
    ]
    for msg in messages:
        print(msg, file=sys.stderr)

    x, y = ks.get_xy(vmin, vmax)
    title = f"{species} {K}-mer histogram"

    ax.bar(x, y, fc="#b2df8a", lw=0)

    if peaks:  # Only works for method 'allpaths'
        t = (ks.min1, ks.max1, ks.min2, ks.max2, ks.min3)
        tcounts = [(x, y) for x, y in ks.counts if x in t]
        if tcounts:
            x, y = zip(*tcounts)
            tcounts = dict(tcounts)
            ax.plot(x, y, "ko", lw=3, mec="k", mfc="w")
            ax.text(ks.max1, tcounts[ks.max1], "SNP peak")
            ax.text(ks.max2, tcounts[ks.max2], "Main peak")

    _, ymax = ax.get_ylim()
    ymax *= 7 / 6
    # Plot the negative binomial fit
    if method == "nbinom":
        plot_nbinom_fit(ax, ks, ymax, method_info)
        messages += [ks.ploidy_message] + ks.copy_messages

    write_messages(ax, messages)

    ax.set_title(markup(title))
    ax.set_xlim((0, vmax))
    ax.set_ylim((0, ymax))
    adjust_spines(ax, ["left", "bottom"], outward=True)
    xlabel, ylabel = "Coverage (X)", "Counts"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    set_human_axis(ax)

    return Genome_size


def histogram(args):
    """
    %prog histogram meryl.histogram species K

    Plot the histogram based on Jellyfish or meryl K-mer distribution, species and N are
    only used to annotate the graphic.
    """
    p = OptionParser(histogram.__doc__)
    p.add_argument(
        "--vmin",
        dest="vmin",
        default=2,
        type=int,
        help="minimum value, inclusive",
    )
    p.add_argument(
        "--vmax",
        dest="vmax",
        default=200,
        type=int,
        help="maximum value, inclusive",
    )
    p.add_argument(
        "--method",
        choices=("nbinom", "allpaths"),
        default="nbinom",
        help="'nbinom' - slow but more accurate for het or polyploid genome; "
        + "'allpaths' - fast and works for homozygous enomes",
    )
    p.add_argument(
        "--maxiter",
        default=100,
        type=int,
        help="Max iterations for optimization. Only used with --method nbinom",
    )
    p.add_argument(
        "--coverage", default=0, type=int, help="Kmer coverage [default: auto]"
    )
    p.add_argument(
        "--nopeaks",
        default=False,
        action="store_true",
        help="Do not annotate K-mer peaks",
    )
    opts, args, iopts = p.set_image_options(args, figsize="7x7")

    if len(args) != 3:
        sys.exit(not p.print_help())

    histfile, species, N = args
    method = opts.method
    vmin, vmax = opts.vmin, opts.vmax
    peaks = not opts.nopeaks and method == "allpaths"
    N = int(N)

    fig = plt.figure(1, (iopts.w, iopts.h))
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))

    Genome_size = draw_ks_histogram(
        ax, histfile, method, opts.coverage, vmin, vmax, species, N, opts.maxiter, peaks
    )

    imagename = histfile.split(".")[0] + "." + iopts.format
    savefig(imagename, dpi=100)

    return Genome_size


if __name__ == "__main__":
    main()
