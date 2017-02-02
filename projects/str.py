#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Related scripts for the HLI-STR (TREDPARSE) paper.
"""

import sys
import pandas as pd

from jcvi.graphics.base import plt, savefig, normalize_axes
from jcvi.apps.base import OptionParser, ActionDispatcher


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
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 0:
        sys.exit(not p.print_help())

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


def compare(args):
    """
    %prog datafile

    Illustrate blablabla...
    """
    p = OptionParser(__doc__)
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
