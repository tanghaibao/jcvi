#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Helper functions for Copy Number Variations (CNV).
"""

import sys
import os.path as op

import numpy as np
import pandas as pd

from collections import defaultdict

from jcvi.utils.aws import sync_from_s3
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, sh


autosomes = ["chr{}".format(x) for x in range(1, 23)]
sexsomes = ["chrX", "chrY"]
allsomes = autosomes + sexsomes


def main():

    actions = (
        ('gcshift', 'correct cib according to GC content'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def load_cib(cibfile, n=1000):
    if not op.exists(cibfile):
        cibfile += ".gz"
    if cibfile.endswith(".gz"):
        sh("pigz -d -k {}".format(cibfile))
        cibfile = cibfile.replace(".gz", "")
    cib = np.fromfile(cibfile, dtype=np.int8) + 128
    rm = pd.rolling_mean(cib, n, min_periods=n / 2)
    a = rm[n - 1::n].copy()
    del cib
    del rm
    return a


def build_gc_array(fastafile="/mnt/ref/hg38.upper.fa",
                   gcdir="gc", n=1000):
    from pyfasta import Fasta
    f = Fasta(fastafile)
    mkdir(gcdir)
    for seqid in allsomes:
        c = np.array(f[seqid])
        gc = (c == 'G') | (c == 'C')  # If base is GC
        rr = ~(c == 'N')              # If base is real
        mgc = pd.rolling_sum(gc, n, min_periods=n / 2)[n - 1::n]
        mrr = pd.rolling_sum(rr, n, min_periods=n / 2)[n - 1::n]
        gc_pct = np.rint(mgc * 100 / mrr)
        gc_pct = np.asarray(gc_pct, dtype=np.uint8)
        arfile = op.join(gcdir, "{}.{}.gc".format(seqid, n))
        gc_pct.tofile(arfile)
        print >> sys.stderr, seqid, gc_pct, arfile


def gcshift(args):
    """
    %prog gcshift \
        s3://hli-bix-us-west-2/kubernetes/wf-root-test/102340_NA12878/lpierce-ccn_gcn-v2/ \
        102340_NA12878

    Download CCN output folder and convert cib to copy number per 1Kb.
    """
    p = OptionParser(gcshift.__doc__)
    p.add_option("--binsize", default=1000, type="int",
                 help="Window size along chromosome")
    p.add_option("--cleanup", default=False, action="store_true",
                 help="Clean up downloaded s3 folder")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    s3dir, sample_key = args
    n = opts.binsize
    gcdir = "gc"
    if not op.exists(gcdir):
        build_gc_array(n=n)

    sync_from_s3(s3dir, target_dir=sample_key)
    # Build GC correction table
    gc_bin = defaultdict(list)
    gc_med = {}
    coverage = []

    for seqid in allsomes:
        gcfile = op.join(gcdir, "{}.{}.gc".format(seqid, n))
        gc = np.fromfile(gcfile, dtype=np.uint8)
        cibfile = op.join(sample_key, "{}.{}.cib".format(sample_key, seqid))
        cib = load_cib(cibfile)
        print >> sys.stderr, seqid, gc.shape[0], cib.shape[0]
        if seqid in autosomes:
            for gci, k in zip(gc, cib):
                gc_bin[gci].append(k)
        coverage.append((seqid, gc, cib))

    for gci, k in gc_bin.items():
        nonzero_k = [x for x in k if x]
        gc_med[gci] = med = np.median(nonzero_k) / 2
        print >> sys.stderr, gci, len(nonzero_k), med

    cndir = sample_key + "-cn"
    mkdir(cndir)
    apply_fun = np.vectorize(gc_med.get)
    # Apply the GC correction over coverage
    for seqid, gc, cib in coverage:
        nitems = cib.shape[0]
        beta = apply_fun(gc[:nitems])
        cn = cib / beta
        cnfile = op.join(cndir, "{}.{}.cn".format(sample_key, seqid))
        cn.tofile(cnfile)

    if opts.cleanup:
        import shutil
        shutil.rmtree(sample_key)


if __name__ == '__main__':
    main()
