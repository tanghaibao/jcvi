#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Processing fastq files
"""

import os.path as op
import sys
import logging

from optparse import OptionParser

from jcvi.utils.iter import pairwise
from jcvi.apps.base import ActionDispatcher, debug, set_grid, sh
debug()

qual_offset = lambda x: 33 if x == "sanger" else 64


class FastqLite (object):
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual

    def __str__(self):
        return "\n".join((self.name, self.seq, "+", self.qual))


class FastqRecord (object):
    def __init__(self, fh, offset=0, key=None):
        self.name = fh.readline().split()
        if not self.name:
            return
        self.name = self.name[0]
        self.seq = fh.readline().rstrip()
        self.l3 = fh.readline().rstrip()
        self.qual = fh.readline().rstrip()
        if offset != 0:
            self.qual = "".join(chr(ord(x) + offset) for x in self.qual)
        self.length = len(self.seq)
        assert self.length == len(self.qual), \
                "length mismatch: seq(%s) and qual(%s)" % (self.seq, self.qual)
        self.id = key(self.name) if key else self.name

    def __str__(self):
        return "\n".join((self.name, self.seq, "+", self.qual))

    def __len__(self):
        return self.length


def iter_fastq(fh, offset=0, key=None):
    logging.debug("reading file `%s`" % fh.name)
    while True:
        rec = FastqRecord(fh, offset=offset, key=key)
        if not rec.name:
            break
        yield rec
    yield None  # guardian


def main():

    actions = (
        ('size', 'total base pairs in the fastq files'),
        ('splitread', 'split appended reads (from JGI)'),
        ('pair', 'pair up two fastq files and combine pairs'),
        ('unpair', 'unpair pairs.fastq files into 1.fastq and 2.fastq'),
        ('pairinplace', 'starting from fragment.fasta, find if ' +\
                "adjacent records can form pairs"),
        ('convert', 'convert between illumina and sanger offset'),
        ('trim', 'trim reads using fastx_trimmer'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def trim(args):
    """
    %prog trim fastqfile

    Wraps `fastx_trimmer` to trim from begin or end of reads.
    """
    p = OptionParser(trim.__doc__)
    set_grid(p)

    p.add_option("-f", dest="first", default=0, type="int",
            help="First base to keep. Default is 1.")
    p.add_option("-l", dest="last", default=0, type="int",
            help="Last base to keep. Default is entire read.")

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    grid = opts.grid

    fastqfile, = args
    base = op.basename(fastqfile).split(".")[0]
    fq = base + ".ntrimmed.fastq"

    cmd = "fastx_trimmer -Q33 "
    if opts.first:
        cmd += "-f {0.first} ".format(opts)
    if opts.last:
        cmd += "-l {0.last} ".format(opts)

    sh(cmd, grid=grid, infile=fastqfile, outfile=fq)


def splitread(args):
    """
    %prog splitread fastqfile

    Split fastqfile into two read fastqfiles, cut in the middle.
    """
    p = OptionParser(splitread.__doc__)
    p.add_option("-n", dest="n", default=76, type="int",
            help="split at N-th base position [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    pairsfastq, = args
    fh = open(pairsfastq)

    assert op.exists(pairsfastq)
    base = op.basename(pairsfastq).split(".")[0]
    fq = base + ".splitted.fastq"
    fw = open(fq, "w")

    it = iter_fastq(fh)
    rec = it.next()
    n = opts.n
    while rec:
        name = rec.name
        seq = rec.seq
        qual = rec.qual

        rec1 = FastqLite(name, seq[:n], qual[:n])
        rec2 = FastqLite(name, seq[n:], qual[n:])

        print >> fw, rec1
        print >> fw, rec2
        rec = it.next()

    logging.debug("reads split into `{0}`".format(fq))
    for f in (fh, fw):
        f.close()


def size(args):
    """
    %prog size fastqfile

    Find the total base pairs in a list of fastq files
    """
    p = OptionParser(size.__doc__)
    opts, args = p.parse_args(args)

    total_size = 0
    total_numrecords = 0
    for f in args:
        fh = open(f)
        for rec in iter_fastq(fh):
            if rec:
                total_numrecords += 1
                total_size += len(rec)

    print >>sys.stderr, "A total %d bases in %s sequences" % (total_size,
            total_numrecords)


def is_low_complexity(seq, cutoff_pct=95):
    # Here we remove seq with one dominant nucleotide
    acgt = defaultdict(int)
    for s in seqs:
        acgt[s] += 1

    for x in acgt.values():
        if x * 100 / len(seq) > cutoff_pct:
            return True

    return False


def convert(args):
    """
    %prog convert in.fastq out.fastq

    illumina fastq quality encoding uses offset 64, and sanger uses 33. This
    script creates a new file with the correct encoding
    """
    p = OptionParser(convert.__doc__)
    p.add_option("-Q", dest="infastq", default="sanger",
            help="input fastq [default: %default]")
    p.add_option("-q", dest="outfastq", default="sanger",
            help="output fastq format [default: %default]")
    p.add_option("-L", dest="remove_lowcomplexity", default=False,
            action="store_true", help="remove poly-A/C/G/T")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    infastq, outfastq = args
    in_offset = qual_offset(opts.infastq)
    out_offset = qual_offset(opts.outfastq)
    offset = out_offset - in_offset

    logging.debug("convert from `%s (%s)` to `%s (%s)`" % (infastq,
        opts.infastq, outfastq, opts.outfastq))

    from jcvi.apps.console import ProgressBar
    totalsize = op.getsize(infastq)
    bar = ProgressBar(maxval=totalsize).start()
    fw = open(outfastq, "w")
    ah = open(infastq)

    low_complexity = 0
    for a in iter_fastq(ah, offset=offset):
        if a is None:
            continue
        if opts.remove_lowcomplexity:
            if is_low_complexity(a.seq):
                low_complexity += 1
                continue

        print >>fw, a

        # update progress
        pos = ah.tell()
        bar.update(pos)

    bar.finish()
    print >>sys.stderr, \
            "A total of %d low complexity region removed." % low_complexity
    sys.stdout.write("\n")


def unpair(args):
    """
    %prog unpair pairs.fastq

    Reverse operation of `pair`:
    The /1 will be placed in 1.fastq, and /2 will be place in 2.fastq.
    """
    p = OptionParser(unpair.__doc__)
    p.add_option("--tag", dest="tag", default=False, action="store_true",
            help="add tag (/1, /2) to the read name")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    pairsfastq, = args
    fh = open(pairsfastq)

    assert op.exists(pairsfastq)
    base = op.basename(pairsfastq).split(".")[0]
    afastq = base + ".1.fastq"
    bfastq = base + ".2.fastq"
    afw = open(afastq, "w")
    bfw = open(bfastq, "w")

    tag = opts.tag

    it = iter_fastq(fh)
    rec = it.next()
    while rec:
        if tag:
            rec.name += "/1"
        print >> afw, rec
        rec = it.next()
        if tag:
            rec.name += "/2"
        print >> bfw, rec
        rec = it.next()

    logging.debug("reads unpaired into `%s` and `%s`" % (afastq, bfastq))
    for f in (fh, afw, bfw):
        f.close()


def pairinplace(args):
    """
    %prog pairinplace bulk.fastq

    Pair up the records in bulk.fastq by comparing the names for adjancent
    records. If they match, print to bulk.pairs.fastq, else print to
    bulk.fragments.fastq.
    """
    p = OptionParser(pairinplace.__doc__)
    p.add_option("-r", dest="rclip", default=0, type="int",
            help="pair ID is derived from rstrip N chars [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    fastqfile, = args
    base = op.basename(fastqfile).split(".")[0]
    frags = base + ".fragments.fastq"
    pairs = base + ".pairs.fastq"

    fragsfw = open(frags, "w")
    pairsfw = open(pairs, "w")

    N = opts.rclip
    if N:
        strip_name = lambda x: x[:-N]
    else:
        strip_name = str

    fh = open(fastqfile)
    fh_iter = iter_fastq(fh, key=strip_name)
    skipflag = False  # controls the iterator skip
    for a, b in pairwise(fh_iter):
        if b is None:  # hit the eof
            break

        if skipflag:
            skipflag = False
            continue

        if a.id == b.id:
            print >> pairsfw, a
            print >> pairsfw, b
            skipflag = True
        else:
            print >> fragsfw, a

    # don't forget the last one, when b is None
    if not skipflag:
        print >> fragsfw, a

    logging.debug("reads paired into `%s` and `%s`" % (pairs, frags))
    for f in (fh, pairsfw, fragsfw):
        f.close()


def pair(args):
    """
    %prog pair 1.fastq 2.fastq

    Pair up the records in 1.fastq and 2.fastq, pairs are indicated by trailing
    "/1" and "/2". If using raw sequences, this is trivial, since we can just
    iterate one at a time for both files; however if two files do not match,
    (e.g. due to trimming), we need a third fastq that provides the order. Two
    output files will be automatically written, one `fragments.fastq` and
    `pairs.fastq`
    """
    p = OptionParser(pair.__doc__)
    p.add_option("-Q", dest="infastq", default="sanger",
            help="input fastq [default: %default]")
    p.add_option("-q", dest="outfastq", default="sanger",
            help="output fastq format [default: %default]")
    p.add_option("-r", dest="ref", default=None,
            help="a reference fastq that provides order")
    p.add_option("-o", dest="outputdir", default=None,
            help="deposit output files into specified directory")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    afastq, bfastq = args

    assert op.exists(afastq) and op.exists(bfastq)
    logging.debug("pair up `%s` and `%s`" % (afastq, bfastq))

    base = op.basename(afastq).split(".")[0]
    frags = base + ".fragments.fastq"
    pairs = base + ".pairs.fastq"

    outputdir = opts.outputdir
    if outputdir:
        frags = op.join(outputdir, frags)
        pairs = op.join(outputdir, pairs)

    logging.debug("output files are `%s` and `%s`" % (frags, pairs))

    in_offset = qual_offset(opts.infastq)
    out_offset = qual_offset(opts.outfastq)
    offset = out_offset - in_offset
    ref = opts.ref

    if ref:
        rh = open(ref)
        totalsize = op.getsize(ref)
    else:
        totalsize = op.getsize(afastq)

    from jcvi.apps.console import ProgressBar
    bar = ProgressBar(maxval=totalsize).start()
    strip_name = lambda x: x.rsplit("/", 1)[0]

    ah = open(afastq)
    ah_iter = iter_fastq(ah, offset=offset, key=strip_name)
    bh = open(bfastq)
    bh_iter = iter_fastq(bh, offset=offset, key=strip_name)

    a = ah_iter.next()
    b = bh_iter.next()

    fragsfw = open(frags, "w")
    pairsfw = open(pairs, "w")

    if ref:
        for r in iter_fastq(rh, offset=0, key=strip_name):
            if not a or not b or not r:
                break

            if a.id == b.id:
                print >>pairsfw, a
                print >>pairsfw, b
                a = ah_iter.next()
                b = bh_iter.next()
            elif a.id == r.id:
                print >>fragsfw, a
                a = ah_iter.next()
            elif b.id == r.id:
                print >>fragsfw, b
                b = bh_iter.next()

            # update progress
            pos = rh.tell()
            bar.update(pos)

        # write all the leftovers to frags file
        while a:
            print >>fragsfw, a
            a = ah_iter.next()

        while b:
            print >>fragsfw, b
            b = bh_iter.next()

    else:  # easy case when afile and bfile records are in order
        # TODO: unimplemented
        pass

    bar.finish()
    sys.stdout.write("\n")


if __name__ == '__main__':
    main()
