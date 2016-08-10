#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Helper functions for Copy Number Variations (CNV).
"""

import sys
import logging
import os.path as op

import numpy as np
import pandas as pd

from collections import defaultdict

from jcvi.utils.aws import sync_from_s3
from jcvi.algorithms.formula import get_kmeans
from jcvi.apps.base import OptionParser, ActionDispatcher, getfilesize, mkdir, sh


autosomes = ["chr{}".format(x) for x in range(1, 23)]
sexsomes = ["chrX", "chrY"]
allsomes = autosomes + sexsomes


def main():

    actions = (
        ('gcshift', 'correct cib according to GC content'),
        ('mergecn', 'compile matrix of GC-corrected copy numbers'),
        # Interact with CCN script
        ('batchccn', 'run CCN script in batch'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def batchccn(args):
    """
    %prog batchccn test.csv

    Run CCN script in batch. Write makefile.
    """
    from jcvi.apps.grid import MakeManager

    p = OptionParser(batchccn.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    csvfile, = args
    mm = MakeManager()
    pf = op.basename(csvfile).split(".")[0]
    mkdir(pf)

    header = open(csvfile).next()
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

    Compile matrix of GC-corrected copy numbers. Place a bunch of folders in the
    csv file. Each folder will be scanned, one chromosomes after another.
    """
    p = OptionParser(mergecn.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    csvfile, = args
    samples = [x.replace("-cn", "").strip().strip("/") for x in open(csvfile)]
    betadir = "beta"
    mkdir(betadir)
    for seqid in allsomes:
        names = [op.join(s + "-cn", "{}.{}.cn".format(s, seqid)) \
                    for s in samples]
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
            logging.debug("K-means with {} c0:{} c1:{}".\
                    format(seqid, zero_med, one_med))
            higher_idx = 1 if one_med > zero_med else 0
            # Use the higher mean coverage componen
            arrays = np.array(arrays)[idx == higher_idx]
        arrays = [[x] for x in arrays]
        ar = np.concatenate(arrays)
        print seqid, ar.shape
        rows, columns = ar.shape
        beta = []
        std = []
        for j in xrange(columns):
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


def load_cib(cibfile, n=1000):
    cibgzfile = cibfile + ".gz"
    if not op.exists(cibfile) or getfilesize(cibfile) < getfilesize(cibgzfile):
        cibfile = cibgzfile
    if cibfile.endswith(".gz"):
        sh("pigz -d -k -f {}".format(cibfile))
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
    cndir = sample_key + "-cn"
    if op.exists(cndir):
        logging.debug("Directory {} exists. Skipped.".format(cndir))
        return

    gcdir = "gc"
    if not op.exists(gcdir):
        build_gc_array(n=n)

    if s3dir.startswith("s3://"):
        sync_from_s3(s3dir, target_dir=sample_key)
    assert op.exists(sample_key), "Directory {} doesn't exist!"\
                    .format(sample_key)

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
