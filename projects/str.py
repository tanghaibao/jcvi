#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Related scripts for the HLI-STR (TREDPARSE) paper.
"""

import os.path as op
import sys
import vcf
import logging
import pandas as pd

from jcvi.graphics.base import normalize_axes, panel_labels, plt, savefig
from jcvi.apps.grid import Parallel
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, iglob


class TREDPARSEvcf(object):

    def __init__(self, vcffile):
        samplekey = op.basename(vcffile).split(".")[0]
        reader = vcf.Reader(open(vcffile, "rb"))
        res = "-1/-1"
        for rec in reader:
            sample = rec.samples[0]
            res = sample["GB"]
            break
        print samplekey, res


def main():

    actions = (
        # Prepare data
        ('mergebam', 'merge sets of BAMs to make diploid'),
        # Compile results
        ('batchlobstr', 'run lobSTR on a list of BAMs'),
        ('compilevcf', 'compile vcf outputs into lists'),
        # Plotting
        ('compare', 'compare callers on fake HD patients'),
        ('evidences', 'plot distribution of evidences'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def mergebam(args):
    """
    %prog mergebam dir1 dir2 homo_outdir
    or
    %prog mergebam dir1 dir2/20.bam het_outdir

    Merge sets of BAMs to make diploid. Two modes:
    - Homozygous mode: pair-up the bams in the two folders and merge
    - Heterozygous mode: pair the bams in first folder with a particular bam
    """
    p = OptionParser(mergebam.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    idir1, idir2, outdir = args
    dir1 = [idir1] if idir1.endswith(".bam") else iglob(idir1, "*.bam")
    dir2 = [idir2] if idir2.endswith(".bam") else iglob(idir2, "*.bam")
    nbams1 = len(dir1)
    nbams2 = len(dir2)
    # Make sure more or the same number of bams in first pile
    if nbams1 < nbams2:
        dir1, dir2 = dir2, dir1
    if nbams1 == nbams2:
        logging.debug("Homozygous mode")
    elif nbams1 > nbams2:
        assert nbams2 == 1, "Second pile must contain a single bam"
        dir2 = [idir2] * nbams1

    assert len(dir1) == len(dir2), "Two piles must contain same number of bams"
    cmd = "samtools merge {} {} {} && samtools index {}"
    cmds = []
    mkdir(outdir)
    for a, b in zip(dir1, dir2):
        ia = op.basename(a).split(".")[0]
        ib = op.basename(b).split(".")[0]
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
    p.add_option("--csv", default="hli.20170316.tred.tsv",
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


def parse_results(datafile, exclude=None):
    fp = open(datafile)
    data = []
    for row in fp:
        truth, call = row.split()
        t = long_allele(truth, exclude=exclude),
        c = long_allele(call, exclude=exclude)
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
    p = OptionParser(__doc__)
    p.add_option('--maxinsert', default=300, type="int",
                 help="Maximum number of repeats")
    opts, args, iopts = p.set_image_options(args, figsize="15x5")

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1,
                                        figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=2)

    # Huntington risk allele
    infected_thr = 40
    ref_thr = 19

    # ax1: Multiple callers at lower range
    df = pd.read_csv("Evaluation.csv")
    truth = df["Truth"]

    ax1.plot(truth, df["Manta"], 'bx-')
    ax1.plot(truth, df["Isaac"], 'yo-')
    ax1.plot(truth, df["GATK"], 'md-')
    ax1.plot(truth, df["lobSTR"], 'c+-')
    ax1.plot(truth, truth, 'k--') # to show diagonal

    bbox = {'facecolor': 'tomato', 'alpha': .2, 'ec': 'w'}
    pad = 2
    ax1.axhline(infected_thr, color='tomato')
    ax1.text(max(truth) - pad, infected_thr + pad, 'Risk threshold',
             bbox=bbox, ha="right")
    ax1.axhline(ref_thr, color='tomato')
    ax1.text(max(truth) - pad, ref_thr - pad, 'Reference repeat count',
             bbox=bbox, ha="right", va="top")

    ax1.set_xlabel(r'Num of CAG repeats inserted ($\mathit{h}$)')
    ax1.set_ylabel('Num of CAG repeats called')
    ax1.set_title(r'Simulated haploid $\mathit{h}$')
    ax1.legend(['Manta', 'Isaac', 'GATK', 'lobSTR', 'Truth'], loc='best')

    max_insert = opts.maxinsert
    # ax2: lobSTR vs TREDPARSE with haploid model
    lobstr_results = parse_results("lobstr_results_homo.txt")
    tredparse_results = parse_results("tredparse_results_homo.txt")
    truth = range(10, max_insert + 1)
    lx, ly = zip(*lobstr_results)
    tx, ty = zip(*tredparse_results)
    lrmsd = compute_rmsd(truth, ly)
    trmsd = compute_rmsd(truth, ty)

    ax2.plot(lx, ly, 'c+-')
    ax2.plot(tx, ty, 'gx-')
    ax2.plot(truth, truth, 'k--')

    ax2.set_xlabel(r'Num of CAG repeats inserted ($\mathit{h}$)')
    ax2.set_ylabel('Num of CAG repeats called')
    ax2.set_title(r'Simulated haploid $\mathit{h}$')
    ax2.legend(['lobSTR (RMSD={:.2f})'.format(lrmsd),
                'TREDPARSE (RMSD={:.2f})'.format(trmsd),
                'Truth'], loc='best')

    pad *= 2
    ax2.axhline(infected_thr, color='tomato')
    ax2.text(max(truth) - pad, infected_thr + pad, 'Risk threshold',
             bbox=bbox, ha="right")
    ax2.set_xlim(10, max_insert)
    ax2.set_ylim(10, max_insert)

    # ax3: lobSTR vs TREDPARSE with haploid model
    lobstr_results = parse_results("lobstr_results_het.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het.txt", exclude=20)
    truth = range(10, max_insert + 1)
    lx, ly = zip(*lobstr_results)
    tx, ty = zip(*tredparse_results)
    lrmsd = compute_rmsd(truth, ly)
    trmsd = compute_rmsd(truth, ty)

    ax3.plot(lx, ly, 'c+-')
    ax3.plot(tx, ty, 'gx-')
    ax3.plot(truth, truth, 'k--')

    ax3.set_xlabel(r'Num of CAG repeats inserted ($\mathit{h}$)')
    ax3.set_ylabel('Num of CAG repeats called')
    ax3.set_title(r'Simulated diploid $\mathit{20/h}$')
    ax3.legend(['lobSTR (RMSD={:.2f})'.format(lrmsd),
                'TREDPARSE (RMSD={:.2f})'.format(trmsd),
                'Truth'], loc='best')
    ax3.axhline(infected_thr, color='tomato')
    ax3.text(max(truth) - pad, infected_thr + pad, 'Risk threshold',
             bbox=bbox, ha="right")
    ax3.set_xlim(10, max_insert)
    ax3.set_ylim(10, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 3., 1 - pad, "B"),
                        (2 / 3., 1 - pad, "C")))
    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
