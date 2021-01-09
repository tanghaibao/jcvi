#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Helper functions for Copy Number Variations (CNV).
"""
from __future__ import print_function

import sys
import logging
import os.path as op

import numpy as np
import numpy.ma as ma
import pandas as pd
import pysam

from collections import Counter, defaultdict
from itertools import groupby
from multiprocessing import Pool
from random import choice
from pybedtools import BedTool, cleanup, set_tempdir

from jcvi.algorithms.formula import get_kmeans
from jcvi.apps.grid import MakeManager
from jcvi.utils.aws import glob_s3, push_to_s3, sync_from_s3
from jcvi.utils.cbook import percentage
from jcvi.apps.base import OptionParser, ActionDispatcher, getfilesize, mkdir, popen, sh


autosomes = ["chr{}".format(x) for x in range(1, 23)]
sexsomes = ["chrX", "chrY"]
allsomes = autosomes + sexsomes
# See: http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/
PAR = [("chrX", 10001, 2781479), ("chrX", 155701383, 156030895)]


class CopyNumberSegment(object):
    def __init__(self, chr, rr, tag, mean_cn, realbins, is_PAR=False):
        self.chr = chr
        self.rr = rr
        self.start = rr[0] * 1000
        self.end = rr[1] * 1000
        self.span = self.end - self.start
        self.tag = tag
        self.mean_cn = mean_cn
        self.realbins = realbins
        self.is_PAR = is_PAR

    def __str__(self):
        mb = self.rr / 1000.0
        coords = "{}:{}-{}Mb".format(self.chr, format_float(mb[0]), format_float(mb[1]))
        if self.is_PAR:
            coords += ":PAR"
        msg = "[{}] {} CN={} bins={}".format(
            self.tag, coords, self.mean_cn, self.realbins
        )
        if self.realbins >= 10000:  # Mark segments longer than 10K bins ~ 10Mb
            msg += "*"
        return msg

    @property
    def bedline(self):
        return "\t".join(
            str(x)
            for x in (self.chr, self.start, self.end, self.tag, self.span, self.mean_cn)
        )


class CopyNumberHMM(object):
    def __init__(
        self, workdir, betadir="beta", mu=0.003, sigma=10, step=0.1, threshold=0.2
    ):
        self.model = self.initialize(mu=mu, sigma=sigma, step=step)
        self.workdir = workdir
        self.betadir = betadir
        if not op.exists(betadir):
            sync_from_s3("s3://hli-mv-data-science/htang/ccn/beta", target_dir=betadir)
        self.mu = mu
        self.sigma = sigma
        self.step = step
        self.threshold = threshold

    def run(self, samplekey, chrs=allsomes):
        if isinstance(chrs, str):
            chrs = [chrs]
        allevents = []
        for chr in chrs:
            X, Z, clen, events = self.run_one(samplekey, chr)
            allevents.extend(events)
        return allevents

    def run_one(self, samplekey, chr):
        cov = np.fromfile(
            "{}/{}-cn/{}.{}.cn".format(self.workdir, samplekey, samplekey, chr)
        )
        beta = np.fromfile("beta/{}.beta".format(chr))
        std = np.fromfile("beta/{}.std".format(chr))
        # Check if the two arrays have different dimensions
        clen, blen = cov.shape[0], beta.shape[0]
        tlen = max(clen, blen)
        if tlen > clen:
            cov = np.array(list(cov) + [np.nan] * (tlen - clen))
        elif tlen > blen:
            beta = np.array(list(beta) + [np.nan] * (tlen - blen))
        clen, blen = cov.shape[0], beta.shape[0]
        assert clen == blen, "cov ({}) and correction ({}) not same dimension".format(
            clen, blen
        )
        normalized = cov / beta
        fixed = normalized.copy()
        fixed[np.where(std > self.threshold)] = np.nan
        X = fixed
        Z = self.predict(X)

        med_cn = np.median(fixed[np.isfinite(fixed)])
        print(chr, med_cn)

        # Annotate segments
        segments = self.annotate_segments(Z)
        events = []
        for mean_cn, rr in segments:
            ss = fixed[rr[0] : rr[1]]
            realbins = np.sum(np.isfinite(ss))
            # Determine whether this is an outlier
            segment = self.tag(chr, mean_cn, rr, med_cn, realbins)
            if segment:
                events.append((mean_cn, rr, segment))
        events.sort(key=lambda x: x[-1].start)

        # Send some debug info to screen
        for mean_cn, rr, segment in events:
            print(segment)

        return X, Z, clen, events

    def tag(self, chr, mean_cn, rr, med_cn, realbins, base=2):
        around_0 = around_value(mean_cn, 0)
        around_1 = around_value(mean_cn, 1)
        around_2 = around_value(mean_cn, 2)
        if realbins <= 1:  # Remove singleton bins
            return
        is_PAR = False
        if chr == "chrX":
            start, end = rr
            is_PAR = end < 5000 or start > 155000
            if med_cn < 1.25:  # Male
                # PAR ~ 2, rest ~ 1
                if is_PAR:
                    base = 2
                    if around_2:
                        return
                else:
                    base = 1
                    if around_1:
                        return
            else:
                # All ~ 2
                if around_2:
                    return
        elif chr == "chrY":
            if med_cn < 0.25:  # Female
                base = 0
                if around_0:
                    return
            else:
                base = 1
                if around_1:
                    return
        else:
            if around_2:
                return
        tag = "DUP" if mean_cn > base else "DEL"
        segment = CopyNumberSegment(chr, rr, tag, mean_cn, realbins, is_PAR=False)
        return segment

    def initialize(self, mu, sigma, step):
        from hmmlearn import hmm

        # Initial population probability
        n = int(10 / step)
        startprob = 1.0 / n * np.ones(n)
        transmat = mu * np.ones((n, n))
        np.fill_diagonal(transmat, 1 - (n - 1) * mu)

        # The means of each component
        means = np.arange(0, step * n, step)
        means.resize((n, 1, 1))
        # The covariance of each component
        covars = sigma * np.ones((n, 1, 1))

        # Build an HMM instance and set parameters
        model = hmm.GaussianHMM(n_components=n, covariance_type="full")

        # Instead of fitting it from the data, we directly set the estimated
        # parameters, the means and covariance of the components
        model.startprob_ = startprob
        model.transmat_ = transmat
        model.means_ = means
        model.covars_ = covars
        return model

    def predict(self, X):
        # Handle missing values
        X = ma.masked_invalid(X)
        mask = X.mask
        dX = ma.compressed(X).reshape(-1, 1)
        dZ = self.model.predict(dX)
        Z = np.array([np.nan for i in range(X.shape[0])])
        Z[~mask] = dZ
        Z = ma.masked_invalid(Z)

        return Z * self.step

    def annotate_segments(self, Z):
        """Report the copy number and start-end segment"""
        # We need a way to go from compressed idices to original indices
        P = Z.copy()
        P[~np.isfinite(P)] = -1
        _, mapping = np.unique(np.cumsum(P >= 0), return_index=True)

        dZ = Z.compressed()
        uniq, idx = np.unique(dZ, return_inverse=True)
        segments = []
        for i, mean_cn in enumerate(uniq):
            if not np.isfinite(mean_cn):
                continue
            for rr in contiguous_regions(idx == i):
                segments.append((mean_cn, mapping[rr]))

        return segments

    def plot(
        self, samplekey, chrs=allsomes, color=None, dx=None, ymax=8, ms=2, alpha=0.7
    ):
        from brewer2mpl import get_map
        import matplotlib.pyplot as plt

        props = dict(boxstyle="round", facecolor="wheat", alpha=0.2)

        if isinstance(chrs, str):
            chrs = [chrs]
        f, axs = plt.subplots(1, len(chrs), sharey=True)
        if not isinstance(axs, np.ndarray):
            axs = np.array([axs])
        plt.tight_layout()
        if color is None:
            color = choice(get_map("Set2", "qualitative", 8).mpl_colors)

        for region, ax in zip(chrs, axs):
            chr, start, end = parse_region(region)
            X, Z, clen, events = self.run_one(samplekey, chr)
            ax.plot(X, ".", label="observations", ms=ms, mfc=color, alpha=alpha)
            ax.plot(Z, "k.", label="hidden", ms=6)

            if start is None and end is None:
                ax.set_xlim(0, clen)
            else:
                ax.set_xlim(start / 1000, end / 1000)

            ax.set_ylim(0, ymax)
            ax.set_xlabel("1Kb bins")
            title = "{} {}".format(samplekey.split("_")[1], chr)
            if dx:
                title += " ({})".format(dx)
            ax.set_title(title)

            # The final calls
            yy = 0.9
            abnormal = [x for x in events if x[-1]]
            if len(abnormal) > 5:
                yinterval = 0.02
                size = 10
            else:
                yinterval = 0.05
                size = 12
            for mean_cn, rr, event in events:
                if mean_cn > ymax:
                    continue
                ax.text(np.mean(rr), mean_cn + 0.2, mean_cn, ha="center", bbox=props)
                if event is None:
                    continue
                ax.text(
                    0.5,
                    yy,
                    str(event).rsplit(" ", 1)[0],
                    color="r",
                    ha="center",
                    transform=ax.transAxes,
                    size=size,
                )
                yy -= yinterval

        axs[0].set_ylabel("Copy number")


def parse_region(region):
    if ":" not in region:
        return region, None, None

    chr, start_end = region.split(":")
    start, end = start_end.split("-")
    return chr, int(start), int(end)


def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    (idx,) = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]  # Edit

    # Reshape the result into two columns
    idx.shape = (-1, 2)
    return idx


def format_float(f):
    s = "{:.3f}".format(f)
    return s.rstrip("0").rstrip(".")


def around_value(s, mu, max_dev=0.25):
    return mu - max_dev < s < mu + max_dev


def main():

    actions = (
        ("cib", "convert bam to cib"),
        ("coverage", "plot coverage along chromosome"),
        ("cn", "correct cib according to GC content"),
        ("mergecn", "compile matrix of GC-corrected copy numbers"),
        ("hmm", "run cnv segmentation"),
        # Gene copy number
        ("exonunion", "collapse overlapping exons within the same gene"),
        ("gcn", "gene copy number based on Canvas results"),
        ("summarycanvas", "count different tags in Canvas vcf"),
        # Interact with CCN script
        ("batchccn", "run CCN script in batch"),
        ("batchcn", "run HMM in batch"),
        ("plot", "plot some chromosomes for visual proof"),
        # Benchmark, training, etc.
        ("sweep", "write a number of commands to sweep parameter space"),
        ("compare", "compare cnv output to ground truths"),
        # Plots
        ("gcdepth", "plot GC content vs depth for genomic bins"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def gcdepth(args):
    """
    %prog gcdepth sample_name tag

    Plot GC content vs depth vs genomnic bins. Inputs are mosdepth output:
    - NA12878_S1.mosdepth.global.dist.txt
    - NA12878_S1.mosdepth.region.dist.txt
    - NA12878_S1.regions.bed.gz
    - NA12878_S1.regions.bed.gz.csi
    - NA12878_S1.regions.gc.bed.gz

    A sample mosdepth.sh script might look like:
    ```
    #!/bin/bash
    LD_LIBRARY_PATH=mosdepth/htslib/ mosdepth/mosdepth $1 \\
        bams/$1.bam -t 4 -c chr1 -n --by 1000

    bedtools nuc -fi GRCh38/WholeGenomeFasta/genome.fa \\
        -bed $1.regions.bed.gz \\
        | pigz -c > $1.regions.gc.bed.gz
    ```
    """
    import hashlib
    from jcvi.algorithms.formula import MAD_interval as confidence_interval
    from jcvi.graphics.base import latex, plt, savefig, set2

    p = OptionParser(gcdepth.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    sample_name, tag = args
    # The tag is used to add to title, also provide a random (hashed) color
    coloridx = int(hashlib.sha256(tag).hexdigest(), 16) % len(set2)
    color = set2[coloridx]

    # mosdepth outputs a table that we can use to plot relationship
    gcbedgz = sample_name + ".regions.gc.bed.gz"
    df = pd.read_csv(gcbedgz, delimiter="\t")
    mf = df.loc[:, ("4_usercol", "6_pct_gc")]
    mf.columns = ["depth", "gc"]

    # We discard any bins that are gaps
    mf = mf[(mf["depth"] > 0.001) | (mf["gc"] > 0.001)]

    # Create GC bins
    gcbins = defaultdict(list)
    for i, row in mf.iterrows():
        gcp = int(round(row["gc"] * 100))
        gcbins[gcp].append(row["depth"])
    gcd = sorted((k * 0.01, confidence_interval(v)) for (k, v) in gcbins.items())
    gcd_x, gcd_y = zip(*gcd)
    m, lo, hi = zip(*gcd_y)

    # Plot
    plt.plot(
        mf["gc"],
        mf["depth"],
        ".",
        color="lightslategray",
        ms=2,
        mec="lightslategray",
        alpha=0.1,
    )
    patch = plt.fill_between(
        gcd_x,
        lo,
        hi,
        facecolor=color,
        alpha=0.25,
        zorder=10,
        linewidth=0.0,
        label="Median +/- MAD band",
    )
    plt.plot(gcd_x, m, "-", color=color, lw=2, zorder=20)

    ax = plt.gca()
    ax.legend(handles=[patch], loc="best")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 100)
    ax.set_title("{} ({})".format(latex(sample_name), tag))
    ax.set_xlabel("GC content")
    ax.set_ylabel("Depth")
    savefig(sample_name + ".gcdepth.png")


def exonunion(args):
    """
    %prog exonunion gencode.v26.annotation.exon.bed

    Collapse overlapping exons within the same gene. File
    `gencode.v26.annotation.exon.bed` can be generated by:

    $ zcat gencode.v26.annotation.gtf.gz |  awk 'OFS="\t" {if ($3=="exon")
    {print $1,$4-1,$5,$10,$12,$14,$16,$7}}' | tr -d '";'
    """
    p = OptionParser(exonunion.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gencodebed,) = args
    beds = BedTool(gencodebed)
    # fields[3] is gene_id; fields[6] is gene_name
    for g, gb in groupby(beds, key=lambda x: x.fields[3]):
        gb = BedTool(gb)
        sys.stdout.write(str(gb.sort().merge(c="4,5,6,7", o=",".join(["first"] * 4))))


def get_gain_loss_summary(vcffile):
    """Extract Canvas:GAIN/LOSS/REF/LOH tags"""
    from cyvcf2 import VCF

    counter = Counter()
    for v in VCF(vcffile):
        tag = v.ID.split(":")[1]
        counter[tag] += 1

    return counter


def summarycanvas(args):
    """
    %prog summarycanvas output.vcf.gz

    Generate tag counts (GAIN/LOSS/REF/LOH) of segments in Canvas output.
    """
    p = OptionParser(summarycanvas.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    for vcffile in args:
        counter = get_gain_loss_summary(vcffile)
        pf = op.basename(vcffile).split(".")[0]
        print(
            pf
            + " "
            + " ".join("{}:{}".format(k, v) for k, v in sorted(counter.items()))
        )


def parse_segments(vcffile):
    """Extract all copy number segments from a CANVAS file

    VCF line looks like:
    chr1    788879  Canvas:GAIN:chr1:788880-821005  N       <CNV>   2       q10
    SVTYPE=CNV;END=821005;CNVLEN=32126      RC:BC:CN:MCC    157:4:3:2
    """
    from cStringIO import StringIO
    from cyvcf2 import VCF

    output = StringIO()
    for v in VCF(vcffile):
        chrom = v.CHROM
        start = v.start
        end = v.INFO.get("END") - 1
        (cn,) = v.format("CN")[0]
        print("\t".join(str(x) for x in (chrom, start, end, cn)), file=output)

    beds = BedTool(output.getvalue(), from_string=True)
    return beds


def counter_mean_and_median(counter):
    """Calculate the mean and median value of a counter"""
    if not counter:
        return np.nan, np.nan

    total = sum(v for k, v in counter.items())
    mid = total / 2
    weighted_sum = 0
    items_seen = 0
    median_found = False
    for k, v in sorted(counter.items()):
        weighted_sum += k * v
        items_seen += v
        if not median_found and items_seen >= mid:
            median = k
            median_found = True
    mean = weighted_sum * 1.0 / total
    return mean, median


def counter_format(counter):
    """Pretty print a counter so that it appears as: "2:200,3:100,4:20" """
    if not counter:
        return "na"

    return ",".join("{}:{}".format(*z) for z in sorted(counter.items()))


def gcn(args):
    """
    %prog gcn gencode.v26.exonunion.bed data/*.vcf.gz

    Compile gene copy njumber based on CANVAS results.
    """
    p = OptionParser(gcn.__doc__)
    p.set_cpus()
    p.set_tmpdir(tmpdir="tmp")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    exonbed = args[0]
    canvasvcfs = args[1:]
    tsvfile = opts.outfile
    tmpdir = opts.tmpdir

    mkdir(tmpdir)
    set_tempdir(tmpdir)

    df = vcf_to_df(canvasvcfs, exonbed, opts.cpus)
    for suffix in (".avgcn", ".medcn"):
        df_to_tsv(df, tsvfile, suffix)


def vcf_to_df_worker(arg):
    """Convert CANVAS vcf to a dict, single thread"""
    canvasvcf, exonbed, i = arg
    logging.debug("Working on job {}: {}".format(i, canvasvcf))
    samplekey = op.basename(canvasvcf).split(".")[0].rsplit("_", 1)[0]
    d = {"SampleKey": samplekey}

    exons = BedTool(exonbed)
    cn = parse_segments(canvasvcf)
    overlaps = exons.intersect(cn, wao=True)
    gcn_store = {}
    for ov in overlaps:
        # Example of ov.fields:
        # [u'chr1', u'11868', u'12227', u'ENSG00000223972.5',
        # u'ENST00000456328.2', u'transcribed_unprocessed_pseudogene',
        # u'DDX11L1', u'.', u'-1', u'-1', u'.', u'0']
        gene_name = "|".join((ov.fields[6], ov.fields[3], ov.fields[5]))
        if gene_name not in gcn_store:
            gcn_store[gene_name] = defaultdict(int)

        cn = ov.fields[-2]
        if cn == ".":
            continue
        cn = int(cn)
        if cn > 10:
            cn = 10
        amt = int(ov.fields[-1])
        gcn_store[gene_name][cn] += amt

    for k, v in sorted(gcn_store.items()):
        v_mean, v_median = counter_mean_and_median(v)
        d[k + ".avgcn"] = v_mean
        d[k + ".medcn"] = v_median
    cleanup()
    return d


def vcf_to_df(canvasvcfs, exonbed, cpus):
    """Compile a number of vcf files into tsv file for easy manipulation"""
    df = pd.DataFrame()
    p = Pool(processes=cpus)
    results = []
    args = [(x, exonbed, i) for (i, x) in enumerate(canvasvcfs)]
    r = p.map_async(vcf_to_df_worker, args, callback=results.append)
    r.wait()

    for res in results:
        df = df.append(res, ignore_index=True)
    return df


def df_to_tsv(df, tsvfile, suffix):
    """Serialize the dataframe as a tsv"""
    tsvfile += suffix
    columns = ["SampleKey"] + sorted(x for x in df.columns if x.endswith(suffix))
    tf = df.reindex_axis(columns, axis="columns")
    tf.sort_values("SampleKey")
    tf.to_csv(tsvfile, sep="\t", index=False, float_format="%.4g", na_rep="na")
    print(
        "TSV output written to `{}` (# samples={})".format(tsvfile, tf.shape[0]),
        file=sys.stderr,
    )


def coverage(args):
    """
    %prog coverage *.coverage

    Plot coverage along chromosome. The coverage file can be generated with:
    $ samtools depth a.bam > a.coverage

    The plot is a simple line plot using matplotlib.
    """
    from jcvi.graphics.base import savefig

    p = OptionParser(coverage.__doc__)
    opts, args, iopts = p.set_image_options(args, format="png")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (covfile,) = args
    df = pd.read_csv(covfile, sep="\t", names=["Ref", "Position", "Depth"])

    xlabel, ylabel = "Position", "Depth"
    df.plot(xlabel, ylabel, color="g")

    image_name = covfile + "." + iopts.format
    savefig(image_name)


def plot(args):
    """
    %prog plot workdir sample chr1,chr2

    Plot some chromosomes for visual proof. Separate multiple chromosomes with
    comma. Must contain folder workdir/sample-cn/.
    """
    from jcvi.graphics.base import savefig

    p = OptionParser(plot.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x7", format="png")

    if len(args) != 3:
        sys.exit(not p.print_help())

    workdir, sample_key, chrs = args
    chrs = chrs.split(",")
    hmm = CopyNumberHMM(workdir=workdir)
    hmm.plot(sample_key, chrs=chrs)

    image_name = sample_key + "_cn." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def sweep(args):
    """
    %prog sweep workdir 102340_NA12878

    Write a number of commands to sweep parameter space.
    """
    p = OptionParser(sweep.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    workdir, sample_key = args
    golden_ratio = (1 + 5 ** 0.5) / 2
    cmd = "python -m jcvi.variation.cnv hmm {} {}".format(workdir, sample_key)
    cmd += " --mu {:.5f} --sigma {:.3f} --threshold {:.3f}"
    mus = [0.00012 * golden_ratio ** x for x in range(10)]
    sigmas = [0.0012 * golden_ratio ** x for x in range(20)]
    thresholds = [0.1 * golden_ratio ** x for x in range(10)]
    print(mus, file=sys.stderr)
    print(sigmas, file=sys.stderr)
    print(thresholds, file=sys.stderr)
    for mu in mus:
        for sigma in sigmas:
            for threshold in thresholds:
                tcmd = cmd.format(mu, sigma, threshold)
                print(tcmd)


def compare_worker(arg):
    cnvoutput, truths = arg
    cmd = "intersectBed -f .5 -F .5"
    cmd += " -a {} -b {} | wc -l".format(cnvoutput, truths)
    nlines = int(popen(cmd, debug=False).read())
    target_lines = len([x for x in open(cnvoutput)])
    truths_lines = len([x for x in open(truths)])
    precision = nlines * 100.0 / target_lines
    recall = nlines * 100.0 / truths_lines
    d = "\t".join(
        str(x)
        for x in (
            cnvoutput,
            truths,
            nlines,
            target_lines,
            truths_lines,
            precision,
            recall,
        )
    )
    return d


def compare(args):
    """
    %prog compare NA12878_array_hg38.bed *.seg

    Compare cnv output to known ground truths.
    """
    p = OptionParser(compare.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    truths = args[0]
    cnvoutputs = args[1:]
    cpus = min(len(cnvoutputs), opts.cpus)
    p = Pool(processes=cpus)
    results = []
    files = [(x, truths) for x in cnvoutputs]
    r = p.map_async(compare_worker, files, callback=results.append)
    r.wait()

    for res in results:
        print("\n".join(res))


def bam_to_cib(arg):
    bamfile, seq, samplekey = arg
    bam = pysam.AlignmentFile(bamfile, "rb")
    name, length = seq["SN"], seq["LN"]
    logging.debug("Computing depth for {} (length={})".format(name, length))
    pileup = bam.pileup(name)
    a = np.ones(length, dtype=np.int8) * -128
    for x in pileup:
        a[x.reference_pos] = min(x.nsegments, 255) - 128

    cibfile = op.join(samplekey, "{}.{}.cib".format(samplekey, name))
    a.tofile(cibfile)
    logging.debug("Depth written to `{}`".format(cibfile))


def cib(args):
    """
    %prog cib bamfile samplekey

    Convert BAM to CIB (a binary storage of int8 per base).
    """
    p = OptionParser(cib.__doc__)
    p.add_option("--prefix", help="Report seqids with this prefix only")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bamfile, samplekey = args
    mkdir(samplekey)
    bam = pysam.AlignmentFile(bamfile, "rb")
    refs = [x for x in bam.header["SQ"]]
    prefix = opts.prefix
    if prefix:
        refs = [x for x in refs if x["SN"].startswith(prefix)]

    task_args = []
    for r in refs:
        task_args.append((bamfile, r, samplekey))
    cpus = min(opts.cpus, len(task_args))
    logging.debug("Use {} cpus".format(cpus))

    p = Pool(processes=cpus)
    for res in p.imap(bam_to_cib, task_args):
        continue


def batchcn(args):
    """
    %prog batchcn workdir samples.csv

    Run CNV segmentation caller in batch mode. Scans a workdir.
    """
    p = OptionParser(batchcn.__doc__)
    p.add_option(
        "--upload",
        default="s3://hli-mv-data-science/htang/ccn",
        help="Upload cn and seg results to s3",
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    workdir, samples = args
    upload = opts.upload
    store = upload + "/{}/*.seg".format(workdir)
    computed = [op.basename(x).split(".")[0] for x in glob_s3(store)]
    computed = set(computed)

    # Generate a bunch of cn commands
    fp = open(samples)
    nskipped = ntotal = 0
    cmd = "python -m jcvi.variation.cnv cn --hmm --cleanup {}".format(workdir)
    for row in fp:
        samplekey, path = row.strip().split(",")
        ntotal += 1
        if samplekey in computed:
            nskipped += 1
            continue
        print(" ".join((cmd, samplekey, path)))

    logging.debug("Skipped: {}".format(percentage(nskipped, ntotal)))


def hmm(args):
    """
    %prog hmm workdir sample_key

    Run CNV segmentation caller. The workdir must contain a subfolder called
    `sample_key-cn` that contains CN for each chromosome. A `beta` directory
    that contains scaler for each bin must also be present in the current
    directory.
    """
    p = OptionParser(hmm.__doc__)
    p.add_option("--mu", default=0.003, type="float", help="Transition probability")
    p.add_option(
        "--sigma",
        default=0.1,
        type="float",
        help="Standard deviation of Gaussian emission distribution",
    )
    p.add_option(
        "--threshold",
        default=1,
        type="float",
        help="Standard deviation must be < this in the baseline population",
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    workdir, sample_key = args
    model = CopyNumberHMM(
        workdir=workdir, mu=opts.mu, sigma=opts.sigma, threshold=opts.threshold
    )
    events = model.run(sample_key)
    params = ".mu-{}.sigma-{}.threshold-{}".format(opts.mu, opts.sigma, opts.threshold)
    hmmfile = op.join(workdir, sample_key + params + ".seg")
    fw = open(hmmfile, "w")
    nevents = 0
    for mean_cn, rr, event in events:
        if event is None:
            continue
        print(" ".join((event.bedline, sample_key)), file=fw)
        nevents += 1
    fw.close()
    logging.debug(
        "A total of {} aberrant events written to `{}`".format(nevents, hmmfile)
    )
    return hmmfile


def batchccn(args):
    """
    %prog batchccn test.csv

    Run CCN script in batch. Write makefile.
    """
    p = OptionParser(batchccn.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (csvfile,) = args
    mm = MakeManager()
    pf = op.basename(csvfile).split(".")[0]
    mkdir(pf)

    header = next(open(csvfile))
    header = None if header.strip().endswith(".bam") else "infer"
    logging.debug("Header={}".format(header))
    df = pd.read_csv(csvfile, header=header)
    cmd = "perl /mnt/software/ccn_gcn_hg38_script/ccn_gcn_hg38.pl"
    cmd += " -n {} -b {}"
    cmd += " -o {} -r hg38".format(pf)
    for i, (sample_key, bam) in df.iterrows():
        cmdi = cmd.format(sample_key, bam)
        outfile = "{}/{}/{}.ccn".format(pf, sample_key, sample_key)
        mm.add(csvfile, outfile, cmdi)
    mm.write()


def mergecn(args):
    """
    %prog mergecn FACE.csv

    Compile matrix of GC-corrected copy numbers. Place a bunch of folders in
    csv file. Each folder will be scanned, one chromosomes after another.
    """
    p = OptionParser(mergecn.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (csvfile,) = args
    samples = [x.replace("-cn", "").strip().strip("/") for x in open(csvfile)]
    betadir = "beta"
    mkdir(betadir)
    for seqid in allsomes:
        names = [
            op.join(s + "-cn", "{}.{}.cn".format(op.basename(s), seqid))
            for s in samples
        ]
        arrays = [np.fromfile(name, dtype=np.float) for name in names]
        shapes = [x.shape[0] for x in arrays]
        med_shape = np.median(shapes)
        arrays = [x for x in arrays if x.shape[0] == med_shape]
        ploidy = 2 if seqid not in ("chrY", "chrM") else 1
        if seqid in sexsomes:
            chr_med = [np.median([x for x in a if x > 0]) for a in arrays]
            chr_med = np.array(chr_med)
            idx = get_kmeans(chr_med, k=2)
            zero_med = np.median(chr_med[idx == 0])
            one_med = np.median(chr_med[idx == 1])
            logging.debug(
                "K-means with {} c0:{} c1:{}".format(seqid, zero_med, one_med)
            )
            higher_idx = 1 if one_med > zero_med else 0
            # Use the higher mean coverage componen
            arrays = np.array(arrays)[idx == higher_idx]
        arrays = [[x] for x in arrays]
        ar = np.concatenate(arrays)
        print(seqid, ar.shape)
        rows, columns = ar.shape
        beta = []
        std = []
        for j in range(columns):
            a = ar[:, j]
            beta.append(np.median(a))
            std.append(np.std(a) / np.mean(a))
        beta = np.array(beta) / ploidy
        betafile = op.join(betadir, "{}.beta".format(seqid))
        beta.tofile(betafile)
        stdfile = op.join(betadir, "{}.std".format(seqid))
        std = np.array(std)
        std.tofile(stdfile)
        logging.debug("Written to `{}`".format(betafile))
        ar.tofile("{}.bin".format(seqid))


def is_matching_gz(origfile, gzfile):
    if not op.exists(origfile):
        return False
    if not op.exists(gzfile):
        return False
    return getfilesize(origfile) == getfilesize(gzfile)


def load_cib(cibfile, n=1000):
    cibgzfile = cibfile + ".gz"
    # When we try unzip if cib not found, or cib does not match cibgz
    if not op.exists(cibfile) or not is_matching_gz(cibfile, cibgzfile):
        if op.exists(cibgzfile):
            cibfile = cibgzfile
    if cibfile.endswith(".gz"):
        sh("pigz -d -k -f {}".format(cibfile))
        cibfile = cibfile.replace(".gz", "")
    if not op.exists(cibfile):
        return

    cib = np.fromfile(cibfile, dtype=np.int8) + 128
    rm = pd.rolling_mean(cib, n, min_periods=n / 2)
    a = rm[n - 1 :: n].copy()
    del cib
    del rm
    return a


def build_gc_array(fastafile="/mnt/ref/hg38.upper.fa", gcdir="gc", n=1000):
    from pyfasta import Fasta

    f = Fasta(fastafile)
    mkdir(gcdir)
    for seqid in allsomes:
        if seqid not in f:
            logging.debug("Seq {} not found. Continue anyway.".format(seqid))
            continue
        c = np.array(f[seqid])
        gc = (c == "G") | (c == "C")  # If base is GC
        rr = ~(c == "N")  # If base is real
        mgc = pd.rolling_sum(gc, n, min_periods=n / 2)[n - 1 :: n]
        mrr = pd.rolling_sum(rr, n, min_periods=n / 2)[n - 1 :: n]
        gc_pct = np.rint(mgc * 100 / mrr)
        gc_pct = np.asarray(gc_pct, dtype=np.uint8)
        arfile = op.join(gcdir, "{}.{}.gc".format(seqid, n))
        gc_pct.tofile(arfile)
        print(seqid, gc_pct, arfile, file=sys.stderr)


def cn(args):
    """
    %prog cn workdir 102340_NA12878 \
        s3://hli-bix-us-west-2/kubernetes/wf-root-test/102340_NA12878/lpierce-ccn_gcn-v2/

    Download CCN output folder and convert cib to copy number per 1Kb.
    """
    p = OptionParser(cn.__doc__)
    p.add_option(
        "--binsize", default=1000, type="int", help="Window size along chromosome"
    )
    p.add_option(
        "--cleanup",
        default=False,
        action="store_true",
        help="Clean up downloaded s3 folder",
    )
    p.add_option(
        "--hmm",
        default=False,
        action="store_true",
        help="Run HMM caller after computing CN",
    )
    p.add_option(
        "--upload",
        default="s3://hli-mv-data-science/htang/ccn",
        help="Upload cn and seg results to s3",
    )
    p.add_option("--rebuildgc", help="Rebuild GC directory rather than pulling from S3")
    opts, args = p.parse_args(args)

    if len(args) == 2:
        workdir, sample_key = args
        s3dir = None
    elif len(args) == 3:
        workdir, sample_key, s3dir = args
    else:
        sys.exit(not p.print_help())

    n = opts.binsize
    rebuildgc = opts.rebuildgc
    mkdir(workdir)
    sampledir = op.join(workdir, sample_key)
    if s3dir:
        sync_from_s3(s3dir, target_dir=sampledir)

    assert op.exists(sampledir), "Directory {} doesn't exist!".format(sampledir)

    cndir = op.join(workdir, sample_key + "-cn")
    if op.exists(cndir):
        logging.debug("Directory {} exists. Skipped.".format(cndir))
        return

    gcdir = "gc"
    if rebuildgc:
        build_gc_array(fastafile=rebuildgc, n=n, gcdir=gcdir)
    if not op.exists(gcdir):
        sync_from_s3("s3://hli-mv-data-science/htang/ccn/gc", target_dir=gcdir)

    # Build GC correction table
    gc_bin = defaultdict(list)
    gc_med = {}
    coverage = []

    for seqid in allsomes:
        gcfile = op.join(gcdir, "{}.{}.gc".format(seqid, n))
        if not op.exists(gcfile):
            logging.error("File {} not found. Continue anyway.".format(gcfile))
            continue
        gc = np.fromfile(gcfile, dtype=np.uint8)
        cibfile = op.join(sampledir, "{}.{}.cib".format(sample_key, seqid))
        cib = load_cib(cibfile)
        print(seqid, gc.shape[0], cib.shape[0], file=sys.stderr)
        if seqid in autosomes:
            for gci, k in zip(gc, cib):
                gc_bin[gci].append(k)
        coverage.append((seqid, gc, cib))

    for gci, k in gc_bin.items():
        nonzero_k = [x for x in k if x]
        gc_med[gci] = med = np.median(nonzero_k) / 2
        print(gci, len(nonzero_k), med, file=sys.stderr)

    mkdir(cndir)
    apply_fun = np.vectorize(gc_med.get)
    # Apply the GC correction over coverage
    for seqid, gc, cib in coverage:
        nitems = cib.shape[0]
        beta = apply_fun(gc[:nitems])
        beta_cn = cib / beta
        cnfile = op.join(cndir, "{}.{}.cn".format(sample_key, seqid))
        beta_cn.tofile(cnfile)

    # Run HMM caller if asked
    segfile = hmm([workdir, sample_key]) if opts.hmm else None

    upload = opts.upload
    if upload:
        push_to_s3(upload, cndir)
        if segfile:
            push_to_s3(upload, segfile)

    if opts.cleanup:
        import shutil

        shutil.rmtree(sampledir)
        shutil.rmtree(cndir)


if __name__ == "__main__":
    main()
