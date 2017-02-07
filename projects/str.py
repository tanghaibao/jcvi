#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Related scripts for the HLI-STR (TREDPARSE) paper.
"""

import sys
import pandas as pd

from jcvi.graphics.base import plt, savefig, normalize_axes
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir


def main():

    actions = (
        ('compare', 'compare callers on fake HD patients'),
        ('evidences', 'plot distribution of evidences'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def evidences(args):
    """
    %prog evidences

    Plot distribution of evidences against two factors:
    - Sample mean coverage
    - Longer allele
    """
    p = OptionParser(__doc__)
    opts, args, iopts = p.set_image_options(args, format="pdf")

    if len(args) != 0:
        sys.exit(not p.print_help())

    format = iopts.format

    # Extract sample coverage first
    df = pd.read_csv("qc-export-MeanCoverage.csv", header=None,
                     names=["Samplekey", "MeanCoverage"], index_col=0)

    # Find coverage for HD
    xf = pd.read_csv("hli.20170201.tred.tsv", sep="\t", index_col=0)
    dp = {}
    tred = "HD"
    for sk, row in xf.iterrows():
        sk = sk.split('_')[1]
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


def compare(args):
    """
    %prog compare Evaluation.csv

    Compare performances of various variant callers on simulated STR datasets.
    """
    p = OptionParser(__doc__)
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]

    # Huntington risk allele
    infected_thr = 40
    ref_thr = 19

    # Read evaluation results
    df = pd.read_csv("Evaluation.csv")
    truth = df["Truth"]

    # Start plotting
    plt.figure(1, (iopts.w, iopts.h))
    plt.plot(truth, df["Manta"], 'bx-')
    plt.plot(truth, df["Isaac"], 'yo-')
    plt.plot(truth, df["GATK"], 'md-')
    plt.plot(truth, df["lobSTR"], 'c+-')
    plt.plot(truth, truth, 'k--') # to show diagonal

    bbox = {'facecolor': 'tomato', 'alpha': .2, 'ec': 'w'}
    pad = 1
    plt.axhline(infected_thr, color='tomato')
    plt.text(min(truth) + pad, infected_thr - pad, 'Risk threshold',
             bbox=bbox, va="top")
    plt.axhline(ref_thr, color='tomato')
    plt.text(max(truth) - pad, ref_thr - pad, 'Reference repeat count',
             bbox=bbox, ha="right", va="top")

    plt.xlabel('Num of CAG repeats inserted')
    plt.ylabel('Num of CAG repeats called')
    plt.title(r'Simulated haploid $\mathit{h}$')
    plt.legend(['Manta', 'Isaac', 'GATK', 'lobSTR', 'Truth'], loc='best')

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
