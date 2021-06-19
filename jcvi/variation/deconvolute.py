#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Deconvolute fastq files according to barcodes.
"""
import os.path as op
import sys
import logging

from collections import namedtuple
from itertools import product, groupby, islice
from multiprocessing import Pool
from more_itertools import flatten

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from jcvi.formats.base import FileMerger, must_open
from jcvi.formats.fastq import FastqPairedIterator
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, glob


def main():

    actions = (
        ("split", "split fastqfile into subsets"),
        ("merge", "consolidate split contents"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


BarcodeLine = namedtuple("BarcodeLine", ["id", "seq"])


def unpack_ambiguous(s):
    """
    List sequences with ambiguous characters in all possibilities.
    """
    sd = [ambiguous_dna_values[x] for x in s]
    return ["".join(x) for x in list(product(*sd))]


def is_barcode_sample(seq, barcode, excludebarcode, trim):
    if seq[:trim] != barcode.seq:
        return False
    hasexclude = any(seq.startswith(x.seq) for x in excludebarcode)
    if hasexclude:
        return False
    return True


def split_barcode_paired(t):

    barcode, excludebarcode, outdir, inputfile = t
    trim = len(barcode.seq)
    outfastq = op.join(outdir, "{0}.{1}.fastq".format(barcode.id, barcode.seq))

    r1, r2 = inputfile
    p1fp, p2fp = FastqPairedIterator(r1, r2)
    fw = open(outfastq, "w")
    while True:
        a = list(islice(p1fp, 4))
        if not a:
            break

        b = list(islice(p2fp, 4))
        title, seq, plus, qual = a
        title, seq, qual = title.strip(), seq.strip(), qual.strip()
        if not is_barcode_sample(seq, barcode, excludebarcode, trim):
            continue

        print("{0}\n{1}\n+\n{2}".format(title, seq[trim:], qual[trim:]), file=fw)
        fw.writelines(b)

    fw.close()


def append_barcode_paired(t):

    barcode, excludebarcode, outdir, inputfile = t
    bs = barcode.seq
    trim = len(bs)
    fake_qual = len(bs) * "#"
    outfastq = op.join(outdir, "{0}.{1}.fastq".format(barcode.id, barcode.seq))

    r1, r2 = inputfile
    p1fp, p2fp = FastqPairedIterator(r1, r2)
    fw = open(outfastq, "w")
    while True:
        a = list(islice(p1fp, 4))
        if not a:
            break

        title, seq, plus, qual = a
        seq = seq.strip()
        if not is_barcode_sample(seq, barcode, excludebarcode, trim):
            continue

        fw.writelines(a)

        title, seq, plus, qual = list(islice(p2fp, 4))
        title, seq, qual = title.strip(), seq.strip(), qual.strip()
        # append barcode
        seq = bs + seq
        qual = fake_qual + qual
        print("{0}\n{1}\n+\n{2}".format(title, seq, qual), file=fw)

    fw.close()


def split_barcode(t):

    barcode, excludebarcode, outdir, inputfile = t
    trim = len(barcode.seq)
    outfastq = op.join(outdir, "{0}.{1}.fastq".format(barcode.id, barcode.seq))

    fp = must_open(inputfile)
    fw = open(outfastq, "w")
    for title, seq, qual in FastqGeneralIterator(fp):
        if not is_barcode_sample(seq, barcode, excludebarcode, trim):
            continue
        print("@{0}\n{1}\n+\n{2}".format(title, seq[trim:], qual[trim:]), file=fw)

    fw.close()


def split(args):
    """
    %prog split barcodefile fastqfile1 ..

    Deconvolute fastq files into subsets of fastq reads, based on the barcodes
    in the barcodefile, which is a two-column file like:
    ID01	AGTCCAG

    Input fastqfiles can be several files. Output files are ID01.fastq,
    ID02.fastq, one file per line in barcodefile.

    When --paired is set, the number of input fastqfiles must be two. Output
    file (the deconvoluted reads) will be in interleaved format.
    """
    p = OptionParser(split.__doc__)
    p.set_outdir(outdir="deconv")
    p.add_option(
        "--nocheckprefix",
        default=False,
        action="store_true",
        help="Don't check shared prefix",
    )
    p.add_option(
        "--paired",
        default=False,
        action="store_true",
        help="Paired-end data",
    )
    p.add_option(
        "--append",
        default=False,
        action="store_true",
        help="Append barcode to 2nd read",
    )
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    barcodefile = args[0]
    fastqfile = args[1:]
    paired = opts.paired
    append = opts.append
    if append:
        assert paired, "--append only works with --paired"

    nfiles = len(fastqfile)

    barcodes = []
    fp = open(barcodefile)
    for row in fp:
        id, seq = row.split()
        for s in unpack_ambiguous(seq):
            barcodes.append(BarcodeLine._make((id, s)))

    nbc = len(barcodes)
    logging.debug("Imported {0} barcodes (ambiguous codes expanded).".format(nbc))
    checkprefix = not opts.nocheckprefix

    if checkprefix:
        # Sanity check of shared prefix
        excludebarcodes = []
        for bc in barcodes:
            exclude = []
            for s in barcodes:
                if bc.id == s.id:
                    continue

                assert bc.seq != s.seq
                if s.seq.startswith(bc.seq) and len(s.seq) > len(bc.seq):
                    logging.error("{0} shares same prefix as {1}.".format(s, bc))
                    exclude.append(s)
            excludebarcodes.append(exclude)
    else:
        excludebarcodes = nbc * [[]]

    outdir = opts.outdir
    mkdir(outdir)

    cpus = opts.cpus
    logging.debug("Create a pool of {0} workers.".format(cpus))
    pool = Pool(cpus)

    if paired:
        assert nfiles == 2, "You asked for --paired, but sent in {0} files".format(
            nfiles
        )
        split_fun = append_barcode_paired if append else split_barcode_paired
        mode = "paired"
    else:
        split_fun = split_barcode
        mode = "single"

    logging.debug("Mode: {0}".format(mode))

    pool.map(
        split_fun, zip(barcodes, excludebarcodes, nbc * [outdir], nbc * [fastqfile])
    )


def merge(args):
    """
    %prog merge folder1 ...

    Consolidate split contents in the folders. The folders can be generated by
    the split() process and several samples may be in separate fastq files. This
    program merges them.
    """
    p = OptionParser(merge.__doc__)
    p.set_outdir(outdir="outdir")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    folders = args
    outdir = opts.outdir
    mkdir(outdir)

    files = flatten(glob("{0}/*.*.fastq".format(x)) for x in folders)
    files = list(files)
    key = lambda x: op.basename(x).split(".")[0]
    files.sort(key=key)
    for id, fns in groupby(files, key=key):
        fns = list(fns)
        outfile = op.join(outdir, "{0}.fastq".format(id))
        FileMerger(fns, outfile=outfile).merge(checkexists=True)


if __name__ == "__main__":
    main()
