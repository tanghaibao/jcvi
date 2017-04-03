#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Related scripts for the HLI-STR (TREDPARSE) paper.
"""

import os.path as op
import os
import sys
import vcf
import logging
import shutil
import numpy as np
import pandas as pd

from pyfaidx import Fasta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcvi.graphics.base import FancyArrow, normalize_axes, panel_labels, plt, savefig
from jcvi.formats.sam import index
from jcvi.variation.str import af_to_counts, read_treds
from jcvi.apps.grid import Parallel
from jcvi.apps.bwa import align
from jcvi.apps.base import sh
from jcvi.assembly.base import wgsim
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, iglob


# Huntington risk allele
infected_thr = 40
ref_thr = 19
SIMULATED_HAPLOID = r'Simulated haploid $\mathit{h}$'
SIMULATED_DIPLOID = r"Simulated diploid $\mathit{20/h}$"


class TREDPARSEvcf(object):

    def __init__(self, vcffile):
        samplekey = op.basename(vcffile).split(".")[0]
        reader = vcf.Reader(open(vcffile, "rb"))
        res = "-1/-1"
        for rec in reader:
            sample = rec.samples[0]
            res = sample["GB"]
            ci = sample["CI"]
            break
        print samplekey, res, ci


def main():

    actions = (
        # Prepare data
        ('simulate', 'simulate bams with varying inserts with dwgsim'),
        ('mergebam', 'merge sets of BAMs to make diploid'),
        # Compile results
        ('batchlobstr', 'run lobSTR on a list of BAMs'),
        ('compilevcf', 'compile vcf outputs into lists'),
        # Plotting
        ('evidences', 'plot distribution of evidences'),
        ('compare', 'compare callers on fake HD patients'),
        ('compare2', 'compare TREDPARSE and lobSTR on fake HD patients'),
        ('compare3', 'compare TREDPARSE on fake HD patients adding evidence'),
        ('compare4', 'compare TREDPARSE on fake HD patients adding coverage'),
        ('allelefreq', 'plot the allele frequencies of some STRs'),
        # Diagram
        ('diagram', 'plot the predictive power of various evidences'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def diagram(args):
    """
    %prog diagram

    Plot the predictive power of various evidences.
    """
    p = OptionParser(diagram.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x4")

    if len(args) != 0:
        sys.exit(not p.print_help())

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Gauge on top, this is log-scale
    lsg = "lightslategray"
    yy = .7
    yinterval = .1
    height = .05
    yp = yy - yinterval - height
    canvas = .95
    xstart = .025
    convert = lambda x: xstart + x * canvas / 1200
    # Symbols
    root.text(.5, .9, r"$L$: Read length, $F$: Flank size, $V$: Pair distance", ha="center")
    root.text(.5, .85, r"ex. $L=150bp, F=9bp, V=500bp$", ha="center")
    root.text(xstart + canvas, yy - height, "STR repeat length", ha="center",
              color=lsg, size=10)

    # Mark the key events
    pad = .02
    arrowlen = canvas * 1.05
    arrowprops = dict(length_includes_head=True, width=.01, fc=lsg, lw=0,
                      head_length=arrowlen * .12, head_width=.04)
    p = FancyArrow(xstart, yy, arrowlen, 0, shape="right", **arrowprops)
    root.add_patch(p)

    ppad = 30
    keyevents = ((       0,               0, -1, r"$0$"),
                 (150 - 18, 150 - 18 - ppad, 0, r"$L - 2F$"),
                  (150 - 9,         150 - 9, 1, r"$L - F$"),
                      (150,      150 + ppad, 2, r"$L$"),
                  (500 - 9,         500 - 9, 3, r"$V - F$"),
             (500 * 2 - 18,    500 * 2 - 18, 2, r"$2(V - F)$"),
                )
    for event, pos, i, label in keyevents:
        _event = convert(event)
        _pos = convert(pos)
        root.plot((_event, _event), (yy - height / 4, yy + height / 4),
                  '-', color='k')
        root.text(_pos, yy + pad, label, rotation=45, va="bottom", size=8)
        if i < 0:
            continue
        ystart = yp - i * yinterval
        root.plot((_event, _event), (ystart, yy - height / 4), ':', color=lsg)

    # Range on bottom. These are simple 4 rectangles, with the range indicating
    # the predictive range.
    CLOSED, OPEN = range(2)
    ranges = ((0,       150 - 18, CLOSED, "Spanning reads"),
              (9,        150 - 9, OPEN,   "Partial reads"),
              (150, 500 * 2 - 18, CLOSED, "Repeat reads"),
              (0,        500 - 9, CLOSED, "Paired-end reads"),
             )
    for start, end, starttag, label in ranges:
        _start = convert(start)
        _end = convert(end)
        data = [[0., 1.], [0., 1.]] if starttag == OPEN else \
               [[1., 0.], [1., 0.]]
        root.imshow(data, interpolation='bicubic', cmap=plt.cm.Greens,
                    extent=[_start, _end, yp, yp + height])
        root.text(_end + pad, yp + height / 2, label, va="center")
        yp -= yinterval

    normalize_axes(root)

    image_name = "diagram." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def plot_allelefreq(ax, df, locus, color='lightslategray'):
    tred = df.ix[locus]
    cnt = af_to_counts(tred["allele_freq"])

    cntx, cnty = zip(*cnt.items())

    motif = tred["motif"]
    cutoff_prerisk = tred["cutoff_prerisk"]
    cutoff_risk = tred["cutoff_risk"]
    npredisease = sum(v for (k, v) in cnt.items() if \
                    cutoff_prerisk <= k < cutoff_risk)
    npatients = sum(v for (k, v) in cnt.items() if k >= cutoff_risk)

    ax.bar(cntx, cnty, fc=color)

    ymin, ymax = ax.get_ylim()
    pad = 1.25
    if cutoff_prerisk < cutoff_risk and npredisease:
        ax.axvline(x=cutoff_prerisk, color="k", lw=2)
        ax.text(cutoff_prerisk + pad, .5 * ymax,
                r"Pre-disease ($\geq${}$\times${}) - {} alleles".\
                format(cutoff_prerisk, motif, npredisease),
                rotation=90, color="k", ha="center", va="center")
    ax.axvline(x=cutoff_risk, color="r", lw=2)
    ax.text(cutoff_risk + pad, .5 * ymax,
            r"Disease ($\geq${}$\times${}) - {} alleles".\
            format(cutoff_risk, motif, npatients),
            rotation=90, color="r", ha="center", va="center")

    x = []  # All allelels
    for k, v in cnt.items():
        x.extend([k] * v)

    ax.set_xlabel("Number of repeat units")
    ax.set_ylabel("Number of alleles")
    ax.set_xlim(0, 50)
    ax.set_title(r"{} ({}) median={:.0f}$\times${}".\
                format(locus, tred["title"], np.median(x), motif))


def allelefreq(args):
    """
    %prog allelefreq HD,DM1,SCA1,SCA17

    Plot the allele frequencies of some STRs.
    """
    p = OptionParser(allelefreq.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 1:
        sys.exit(not p.print_help())

    loci, = args
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2,
                                                 figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=4)
    treds, df = read_treds()
    df = df.set_index(["abbreviation"])

    for ax, locus in zip((ax1, ax2, ax3, ax4), loci.split(",")):
        plot_allelefreq(ax, df, locus)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B"),
                        (pad / 2, 1 / 2. , "C"), (1 / 2., 1 / 2. , "D")))
    normalize_axes(root)

    image_name = "allelefreq." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def make_fasta(seq, fastafile, id):
    rec = SeqRecord(Seq(seq), description="", id=id)
    fw = open(fastafile, "w")
    SeqIO.write([rec], fw, "fasta")
    fw.close()


def add_simulate_options(p):
    p.add_option("--readlen", default=150, type="int",
                 help="Length of the read")
    p.add_option("--distance", default=500, type="int",
                 help="Outer distance between the two ends")
    p.set_depth(depth=20)


def simulate(args):
    """
    %prog simulate run_dir 1 300

    Simulate BAMs with varying inserts with dwgsim. The above command will
    simulate between 1 to 300 CAGs in the HD region, in a directory called
    `run_dir`.
    """
    p = OptionParser(simulate.__doc__)
    p.add_option("--ref", default="/Users/htang/projects/ref/hg38.upper.fa",
                 help="Reference genome sequence")
    add_simulate_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    rundir, startunits, endunits = args
    startunits, endunits = int(startunits), int(endunits)
    basecwd = os.getcwd()
    mkdir(rundir)
    os.chdir(rundir)
    cwd = os.getcwd()

    # Huntington region
    pad_left, pad_right = 1000, 10000
    chr, start, end = 'chr4', 3074877, 3074933
    fasta = Fasta(opts.ref)
    seq_left = fasta[chr][start - pad_left:start - 1]
    seq_right = fasta[chr][end: end + pad_right]
    motif = 'CAG'
    reffastafile = "ref.fasta"
    seq = str(fasta[chr][start - pad_left: end + pad_right])
    make_fasta(seq, reffastafile, id=chr.upper())

    # Write fake sequence
    for units in range(startunits, endunits + 1):
        pf = str(units)
        mkdir(pf)
        os.chdir(pf)
        seq = str(seq_left) + motif * units + str(seq_right)
        fastafile = pf + ".fasta"
        make_fasta(seq, fastafile, id=chr.upper())

        # Simulate reads on it
        wgsim([fastafile, "--depth={}".format(opts.depth),
                          "--readlen={}".format(opts.readlen),
                          "--distance={}".format(opts.distance),
                          "--outfile={}".format(pf)])

        read1 = pf + ".bwa.read1.fastq"
        read2 = pf + ".bwa.read2.fastq"
        samfile, _ = align(["../{}".format(reffastafile), read1, read2])
        indexed_samfile = index([samfile])

        sh("mv {} ../{}.bam".format(indexed_samfile, pf))
        sh("mv {}.bai ../{}.bam.bai".format(indexed_samfile, pf))

        os.chdir(cwd)
        shutil.rmtree(pf)

    os.chdir(basecwd)


def mergebam(args):
    """
    %prog mergebam dir1 homo_outdir
    or
    %prog mergebam dir1 dir2/20.bam het_outdir

    Merge sets of BAMs to make diploid. Two modes:
    - Homozygous mode: pair-up the bams in the two folders and merge
    - Heterozygous mode: pair the bams in first folder with a particular bam
    """
    p = OptionParser(mergebam.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) not in (2, 3):
        sys.exit(not p.print_help())

    if len(args) == 2:
        idir1, outdir = args
        dir1 = [idir1] if idir1.endswith(".bam") else iglob(idir1, "*.bam")
        logging.debug("Homozygous mode")
        dir2 = [""] * len(dir1)
    elif len(args) == 3:
        idir1, idir2, outdir = args
        dir1 = [idir1] if idir1.endswith(".bam") else iglob(idir1, "*.bam")
        dir2 = [idir2] if idir2.endswith(".bam") else iglob(idir2, "*.bam")
        assert len(dir2) == 1, "Second pile must contain a single bam"
        dir2 = [idir2] * len(dir1)

    assert len(dir1) == len(dir2), "Two piles must contain same number of bams"
    cmd = "samtools merge {} {} {} && samtools index {}"
    cmds = []
    mkdir(outdir)
    for a, b in zip(dir1, dir2):
        ia = op.basename(a).split(".")[0]
        ib = op.basename(b).split(".")[0] if b else ia
        outfile = op.join(outdir, "{}_{}.bam".format(ia, ib))
        cmds.append(cmd.format(outfile, a, b, outfile))

    p = Parallel(cmds, cpus=opts.cpus)
    p.run()


def batchlobstr(args):
    """
    %prog batchlobstr bamlist

    Run lobSTR on a list of BAMs. The corresponding batch command for TREDPARSE:
    $ tred.py --toy bamlist --haploid CHR4 --workdir tredparse_results
    """
    p = OptionParser(batchlobstr.__doc__)
    p.add_option("--haploid", default="chrY,chrM",
                 help="Use haploid model for these chromosomes")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bamlist, = args
    cmd = "python -m jcvi.variation.str lobstr TOY"
    cmd += " --input_bam_path {}"
    cmd += " --haploid {}".format(opts.haploid)
    cmds = [cmd.format(x.strip()) for x in open(bamlist).readlines()]
    p = Parallel(cmds, cpus=opts.cpus)
    p.run()


def compilevcf(args):
    """
    %prog compilevcf dir

    Compile vcf outputs into lists.
    """
    from jcvi.variation.str import LobSTRvcf

    p = OptionParser(compilevcf.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    folder, = args
    vcf_files = iglob(folder, "*.vcf,*.vcf.gz")
    for vcf_file in vcf_files:
        try:
            p = LobSTRvcf(columnidsfile=None)
            p.parse(vcf_file, filtered=False)
            res = p.items()
            if res:
                k, v = res[0]
                res = v.replace(',', '/')
            else:
                res = "-1/-1"
            num = op.basename(vcf_file).split(".")[0]
            print num, res
        except (TypeError, AttributeError) as e:
            p = TREDPARSEvcf(vcf_file)
            continue


def evidences(args):
    """
    %prog evidences

    Plot distribution of evidences against two factors:
    - Sample mean coverage
    - Longer allele
    """
    p = OptionParser(evidences.__doc__)
    p.add_option("--csv", default="hli.20170328.tred.tsv",
                 help="TRED csv output to plot")
    opts, args, iopts = p.set_image_options(args, format="pdf")

    if len(args) != 0:
        sys.exit(not p.print_help())

    format = iopts.format

    # Extract sample coverage first
    df = pd.read_csv("qc-export-MeanCoverage.csv", header=None,
                     names=["Samplekey", "MeanCoverage"], index_col=0)

    # Find coverage for HD
    xf = pd.read_csv(opts.csv, sep="\t", index_col=0)
    dp = {}
    tred = "HD"
    for sk, row in xf.iterrows():
        sk = str(sk)
        a1 = row[tred + ".1"]
        a2 = row[tred + ".2"]
        fdp = row[tred + ".FDP"]
        pdp = row[tred + ".PDP"]
        pedp = row[tred + ".PEDP"]
        dp[sk] = (a1, a2, fdp, pdp, pedp)

    # Build a consolidated dataframe
    ef = pd.DataFrame.from_dict(dp, orient="index")
    ef.columns = [tred + ".1", tred + ".2", tred + ".FDP",
                  tred + ".PDP", tred + ".PEDP"]
    ef.index.name = "SampleKey"
    mf = df.merge(ef, how="right", left_index=True, right_index=True)

    # Plot a bunch of figures
    outdir = "output"
    mkdir(outdir)
    xlim = ylim = (0, 100)
    draw_jointplot(outdir + "/A", "MeanCoverage", "HD.FDP",
                   data=mf, xlim=xlim, ylim=ylim, format=format)
    draw_jointplot(outdir + "/B", "MeanCoverage", "HD.PDP",
                   data=mf, color='g', xlim=xlim, ylim=ylim, format=format)
    draw_jointplot(outdir + "/C", "MeanCoverage", "HD.PEDP",
                   data=mf, color='m', xlim=xlim, ylim=ylim, format=format)

    xlim = (0, 50)
    draw_jointplot(outdir + "/D", "HD.2", "HD.FDP",
                   data=mf, xlim=xlim, ylim=ylim, format=format)
    draw_jointplot(outdir + "/E", "HD.2", "HD.PDP",
                   data=mf, color='g', xlim=xlim, ylim=ylim, format=format)
    draw_jointplot(outdir + "/F", "HD.2", "HD.PEDP",
                   data=mf, color='m', xlim=xlim, ylim=ylim, format=format)


def draw_jointplot(figname, x, y, data=None, kind="reg", color=None,
                   xlim=None, ylim=None, format="pdf"):
    """
    Wraps around sns.jointplot
    """
    import seaborn as sns
    sns.set_context('talk')
    plt.clf()

    register = {"MeanCoverage": "Sample Mean Coverage",
                "HD.FDP": "Depth of full spanning reads",
                "HD.PDP": "Depth of partial spanning reads",
                "HD.PEDP": "Depth of paired-end reads",
                "HD.2": "Repeat size of the longer allele"}

    g = sns.jointplot(x, y, data=data, kind=kind, color=color,
                      xlim=xlim, ylim=ylim)
    g.ax_joint.set_xlabel(register.get(x, x))
    g.ax_joint.set_ylabel(register.get(y, y))
    savefig(figname + "." + format, cleanup=False)


def long_allele(s, default=19, exclude=None):
    if '_' in s:
        a, b = s.split('_')
    elif '/' in s:
        a, b = s.split('/')
    else:
        raise Exception, "Don't know how to split string {}".format(s)

    res = [int(a), int(b)]
    if exclude and exclude in res:
        res.remove(exclude)
    res = max(res)
    return default if res < 0 else res


def get_lo_hi_from_CI(s, exclude=None):
    """
    Parse the confidence interval from CI.

    >>> get_lo_hi_from_CI("20-20/40-60")
    (40, 60)
    """
    a, b = s.split("|")
    ai, aj = a.split("-")
    bi, bj = b.split("-")

    los = [int(ai), int(bi)]
    his = [int(aj), int(bj)]
    if exclude and exclude in los:
        los.remove(exclude)
    if exclude and exclude in his:
        his.remove(exclude)
    return max(los), max(his)


def parse_results(datafile, exclude=None):
    fp = open(datafile)
    data = []
    for row in fp:
        atoms = row.split()
        truth, call = atoms[:2]
        t = long_allele(truth, exclude=exclude)
        c = long_allele(call, exclude=exclude)
        if len(atoms) == 3:
            ci = atoms[2]
            lo, hi = get_lo_hi_from_CI(ci, exclude=exclude)
            if lo > c:
                lo = c
            if hi < c:
                hi = c
            data.append((t, c, lo, hi))
        else:
            data.append((t, c))
    return data


def compute_rmsd(truth, a):
    if len(a) > len(truth):
        a = a[: len(truth)]
    return (sum((i - j) ** 2 for (i, j) in zip(truth, a)) / len(truth)) ** .5


def compare(args):
    """
    %prog compare Evaluation.csv

    Compare performances of various variant callers on simulated STR datasets.
    """
    p = OptionParser(compare.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2,
                                                 figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=3)

    bbox = {'facecolor': 'tomato', 'alpha': .2, 'ec': 'w'}
    pad = 2

    # Read benchmark data
    df = pd.read_csv("Evaluation.csv")
    truth = df["Truth"]
    axes = (ax1, ax2, ax3, ax4)
    progs = ("Manta", "Isaac", "GATK", "lobSTR")
    markers = ("bx-", "yo-", "md-", "c+-")

    for ax, prog, marker in zip(axes, progs, markers):
        ax.plot(truth, df[prog], marker)
        ax.plot(truth, truth, 'k--') # to show diagonal
        ax.axhline(infected_thr, color='tomato')
        ax.text(max(truth) - pad, infected_thr + pad, 'Risk threshold',
                bbox=bbox, ha="right")
        ax.axhline(ref_thr, color='tomato')
        ax.text(max(truth) - pad, ref_thr - pad, 'Reference repeat count',
                bbox=bbox, ha="right", va="top")
        ax.set_title(SIMULATED_HAPLOID)
        ax.set_xlabel(r'Num of CAG repeats inserted ($\mathit{h}$)')
        ax.set_ylabel('Num of CAG repeats called')
        ax.legend([prog, 'Truth'], loc='best')

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B"),
                        (pad / 2, 1 / 2.,  "C"), (1 / 2., 1 / 2. , "D")))
    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def plot_compare(ax, title, tredparse_results, lobstr_results, pad=8, ms=3,
                 max_insert=300, depth=20, readlen=150, distance=500, color='g',
                 risk=True):
    truth = range(1, max_insert + 1)
    tx, ty, tl, th = zip(*tredparse_results)
    trmsd = compute_rmsd(truth, ty)
    if lobstr_results:
        lx, ly = zip(*lobstr_results)
        lrmsd = compute_rmsd(truth, ly)

    if lobstr_results:
        ax.plot(lx, ly, 'c+-', ms=ms, label='lobSTR (RMSD={:.2f})'.format(lrmsd))
    ax.plot(tx, ty, '.-', color=color, ms=ms, label='TREDPARSE (RMSD={:.2f})'.format(trmsd))
    ax.plot(truth, truth, 'k--', label='Truth')
    ax.fill_between(tx, tl, th, facecolor=color, alpha=.25,
                     label='TREDPARSE 95$\%$ CI')

    ax.set_xlabel(r'Num of CAG repeats inserted ($\mathit{h}$)')
    ax.set_ylabel('Num of CAG repeats called')
    ax.set_title(title)
    ax.legend(loc='best')

    bbox = {'facecolor': 'tomato', 'alpha': .2, 'ec': 'w'}
    if risk:
        ax.axhline(infected_thr, color='tomato')
        ax.text(max(truth) - pad, infected_thr + pad,
                 'Risk cutoff={}'.format(infected_thr) +
                 r'$\times$CAGs',
                 bbox=bbox, ha="right")
    else:
        readlength, pairdistance = 150 / 3, 500 / 3
        ax.axhline(readlength, color='tomato')
        ax.text(max(truth) - pad, readlength + pad,
                 'Read Length ($L$)',
                 bbox=bbox, ha="right")
        ax.axhline(pairdistance, color='tomato')
        ax.text(max(truth) - pad, pairdistance + pad,
                 'Paired-end distance($V$)',
                 bbox=bbox, ha="right")


def compare2(args):
    """
    %prog compare2

    Compare performances of various variant callers on simulated STR datasets.
    """
    p = OptionParser(compare2.__doc__)
    p.add_option('--maxinsert', default=300, type="int",
                 help="Maximum number of repeats")
    add_simulate_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="10x5")

    if len(args) != 0:
        sys.exit(not p.print_help())

    depth = opts.depth
    readlen = opts.readlen
    distance = opts.distance
    max_insert = opts.maxinsert
    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1,
                                   figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=2)

    # ax1: lobSTR vs TREDPARSE with haploid model
    lobstr_results = parse_results("lobstr_results_homo.txt")
    tredparse_results = parse_results("tredparse_results_homo.txt")
    title = SIMULATED_HAPLOID + \
            r" ($D=%s\times, L=%dbp, V=%dbp$)" % (depth, readlen, distance)
    plot_compare(ax1, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    # ax2: lobSTR vs TREDPARSE with diploid model
    lobstr_results = parse_results("lobstr_results_het.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het.txt", exclude=20)
    title = SIMULATED_DIPLOID + \
            r" ($D=%s\times, L=%dbp, V=%dbp$)" % (depth, readlen, distance)
    plot_compare(ax2, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    for ax in (ax1, ax2):
        ax.set_xlim(0, max_insert)
        ax.set_ylim(0, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B")))
    normalize_axes(root)

    image_name = "tredparse." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def compare3(args):
    """
    %prog compare3

    Compare performances of various variant callers on simulated STR datasets.
    This compares the power of various evidence types.
    """
    p = OptionParser(compare3.__doc__)
    p.add_option('--maxinsert', default=300, type="int",
                 help="Maximum number of repeats")
    add_simulate_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 0:
        sys.exit(not p.print_help())

    max_insert = opts.maxinsert
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2,
                                   figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=3)

    color = "lightslategray"
    # ax1: Spanning
    tredparse_results = parse_results("tredparse_results_het-spanning.txt")
    title = SIMULATED_DIPLOID + "( Sub-model 1: Spanning reads)"
    plot_compare(ax1, title, tredparse_results, None, color=color,
                 max_insert=max_insert, risk=False)

    # ax2: Partial
    tredparse_results = parse_results("tredparse_results_het-partial.txt", exclude=20)
    title = SIMULATED_DIPLOID + " (Sub-model 2: Partial reads)"
    plot_compare(ax2, title, tredparse_results, None, color=color,
                 max_insert=max_insert, risk=False)

    # ax3: Repeat
    tredparse_results = parse_results("tredparse_results_het-repeat.txt", exclude=20)
    # HACK (repeat reads won't work under 50)
    tredparse_results = [x for x in tredparse_results if x[0] > 50]
    title = SIMULATED_DIPLOID + " (Sub-model 3: Repeat-only reads)"
    plot_compare(ax3, title, tredparse_results, None, color=color,
                 max_insert=max_insert, risk=False)

    # ax4: Pair
    tredparse_results = parse_results("tredparse_results_het-pair.txt", exclude=20)
    title = SIMULATED_DIPLOID + " (Sub-model 4: Paired-end reads)"
    plot_compare(ax4, title, tredparse_results, None, color=color,
                 max_insert=max_insert, risk=False)

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlim(0, max_insert)
        ax.set_ylim(0, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B"),
                        (pad / 2, 1 / 2. , "C"), (1 / 2., 1 / 2. , "D")))
    normalize_axes(root)

    image_name = "tredparse." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def compare4(args):
    """
    %prog compare4

    Compare performances of various variant callers on simulated STR datasets.
    Adds coverage comparisons as panel C and D.
    """
    p = OptionParser(compare4.__doc__)
    p.add_option('--maxinsert', default=300, type="int",
                 help="Maximum number of repeats")
    add_simulate_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 0:
        sys.exit(not p.print_help())

    depth = opts.depth
    max_insert = opts.maxinsert
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2,
                                   figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=3)

    # ax1: lobSTR vs TREDPARSE with haploid model
    lobstr_results = parse_results("lobstr_results_homo-20x-150bp-500bp.txt")
    tredparse_results = parse_results("tredparse_results_homo-20x-150bp-500bp.txt")
    title = SIMULATED_HAPLOID + r" ($Depth=%s\times)" % depth
    plot_compare(ax1, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    # ax2: lobSTR vs TREDPARSE with diploid model (depth=20x)
    lobstr_results = parse_results("lobstr_results_het-20x-150bp-500bp.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het-20x-150bp-500bp.txt", exclude=20)
    title = SIMULATED_DIPLOID + r" ($Depth=%s\times$)" % depth
    plot_compare(ax2, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    # ax3: lobSTR vs TREDPARSE with diploid model (depth=5x)
    lobstr_results = parse_results("lobstr_results_het-5x-150bp-500bp.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het-5x-150bp-500bp.txt", exclude=20)
    title = SIMULATED_DIPLOID + r" ($Depth=%s\times$)" % 5
    plot_compare(ax3, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    # ax4: lobSTR vs TREDPARSE with diploid model (depth=80x)
    lobstr_results = parse_results("lobstr_results_het-80x-150bp-500bp.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het-80x-150bp-500bp.txt", exclude=20)
    title = SIMULATED_DIPLOID + r" ($Depth=%s\times$)" % 80
    plot_compare(ax4, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlim(0, max_insert)
        ax.set_ylim(0, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B"),
                        (pad / 2, 1 / 2. , "C"), (1 / 2., 1 / 2. , "D")))
    normalize_axes(root)

    image_name = "tredparse." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
