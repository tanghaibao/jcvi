#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Processing fastq files
"""

import os.path as op
import sys
import logging

from optparse import OptionParser

from jcvi.formats.fasta import must_open
from jcvi.apps.base import getfilesize
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

    @property
    def quality(self):
        return [ord(x) for x in self.qual]


def iter_fastq(filename, offset=0, key=None):
    if isinstance(filename, str):
        logging.debug("Read file `{0}`".format(filename))
        fh = must_open(filename)
    else:
        fh = filename

    while True:
        rec = FastqRecord(fh, offset=offset, key=key)
        if not rec.name:
            break
        yield rec
    yield None  # sentinel


def main():

    actions = (
        ('size', 'total base pairs in the fastq files'),
        ('shuffle', 'shuffle paired reads into the same file interleaved'),
        ('split', 'split paired reads into two files'),
        ('splitread', 'split appended reads (from JGI)'),
        ('pair', 'pair up two fastq files and combine pairs'),
        ('unpair', 'unpair pairs.fastq files into 1.fastq and 2.fastq'),
        ('pairinplace', 'starting from fragment.fasta, find if ' +\
                "adjacent records can form pairs"),
        ('convert', 'convert between illumina and sanger offset'),
        ('trim', 'trim reads using fastx_trimmer'),
        ('some', 'select a subset of fastq reads'),
        ('guessoffset', 'guess the quality offset of the fastq records'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def checkShuffleSizes(p1, p2, pairsfastq):
    from jcvi.apps.base import getfilesize

    pairssize = getfilesize(pairsfastq)
    p1size = getfilesize(p1)
    p2size = getfilesize(p2)
    assert pairssize == p1size + p2size, \
          "The sizes do not add up: {0} + {1} != {2}".\
          format(p1size, p2size, pairssize)


def shuffle(args):
    """
    %prog shuffle p1.fastq p2.fastq pairs.fastq

    Shuffle pairs into interleaved format, using `shuffleSequences_fastq.pl`.
    """
    p = OptionParser(shuffle.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    p1, p2, pairsfastq = args
    cmd = "shuffleSequences_fastq.pl {0} {1} {2}".format(*args)
    sh(cmd)

    checkShuffleSizes(p1, p2, pairsfastq)


def split(args):
    """
    %prog split pairs.fastq

    Split shuffled pairs into `.1.fastq` and `.2.fastq`, using `sed`.
    <http://seqanswers.com/forums/showthread.php?t=13776>
    """
    from jcvi.apps.grid import Jobs

    p = OptionParser(split.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pairsfastq, = args
    pf = pairsfastq.rsplit(".", 1)[0]
    p1 = pf + ".1.fastq"
    p2 = pf + ".2.fastq"
    p1cmd = "sed -ne '1~8{{N;N;N;p}}' {0} > {1}".format(pairsfastq, p1)
    p2cmd = "sed -ne '5~8{{N;N;N;p}}' {0} > {1}".format(pairsfastq, p2)

    args = [(p1cmd, ), (p2cmd, )]
    m = Jobs(target=sh, args=args)
    m.run()

    checkShuffleSizes(p1, p2, pairsfastq)


def guessoffset(args):
    """
    %prog guessoffset fastqfile

    Guess the quality offset of the fastqfile, whether 33 or 64.
    See encoding schemes: <http://en.wikipedia.org/wiki/FASTQ_format>

      SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS...............................
      ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
      LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL...............................
      !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh
      |                         |    |        |                              |
     33                        59   64       73                            104

     S - Sanger        Phred+33,  raw reads typically (0, 40)
     X - Solexa        Solexa+64, raw reads typically (-5, 40)
     I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
     J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
     L - Illumina 1.8+ Phred+33,  raw reads typically (0, 40)
        with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
    """
    p = OptionParser(guessoffset.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    ai = iter_fastq(fastqfile)
    rec = ai.next()
    offset = 64
    while rec:
        quality = rec.quality
        lowcounts = len([x for x in quality if x < 59])
        highcounts = len([x for x in quality if x > 74])
        diff = highcounts - lowcounts
        if diff > 10:
            break
        elif diff < -10:
            offset = 33
            break
        rec = ai.next()

    if offset == 33:
        print >> sys.stderr, "Sanger encoding (offset=33)"
    elif offset == 64:
        print >> sys.stderr, "Illumina encoding (offset=64)"

    return offset


def some(args):
    """
    %prog some idsfile afastq bfastq

    Select a subset of the reads with ids present in the idsfile.
    """
    p = OptionParser(some.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    idsfile, afastq, bfastq = args
    fp = open(idsfile)

    ai = iter_fastq(open(afastq))
    arec = ai.next()
    bi = iter_fastq(open(bfastq))
    brec = bi.next()

    for row in fp:
        name = row.strip()
        while arec:
            if arec.name[1:] == name:
                print arec
                print brec
                break
            else:
                arec = ai.next()
                brec = bi.next()


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
        sys.exit(not p.print_help())

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
        sys.exit(not p.print_help())

    pairsfastq, = args

    assert op.exists(pairsfastq)
    base = op.basename(pairsfastq).split(".")[0]
    fq = base + ".splitted.fastq"
    fw = open(fq, "w")

    it = iter_fastq(pairsfastq)
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
        for rec in iter_fastq(f):
            if rec:
                total_numrecords += 1
                total_size += len(rec)

    print >>sys.stderr, "A total %d bases in %s sequences" % (total_size,
            total_numrecords)


def convert(args):
    """
    %prog convert in.fastq out.fastq

    illumina fastq quality encoding uses offset 64, and sanger uses 33. This
    script creates a new file with the correct encoding
    """
    p = OptionParser(convert.__doc__)
    p.add_option("-Q", dest="infastq", default="illumina",
            help="input fastq [default: %default]")
    p.add_option("-q", dest="outfastq", default="sanger",
            help="output fastq format [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    infastq, outfastq = args

    from jcvi.apps.command import EMBOSSPATH

    seqret = EMBOSSPATH("seqret")
    if infastq.endswith(".gz"):
        cmd = "zcat {0} | ".format(infastq)
        cmd += seqret + " fastq-{0}::stdin fastq-{1}::stdout".\
                format(opts.infastq, opts.outfastq)
        cmd += " | gzip > {0}".format(outfastq)
    else:
        cmd = seqret + " fastq-{0}::{1} fastq-{2}::{3}".\
                format(opts.infastq, infastq, opts.outfastq, outfastq)

    sh(cmd)

    return outfastq


def unpair(args):
    """
    %prog unpair *.pairs.fastq unpaired

    Reverse operation of `pair`:
    /1 will be placed in unpaired.1.fastq,
    /2 will be placed in unpaired.2.fastq.
    """
    p = OptionParser(unpair.__doc__)
    p.add_option("--tag", dest="tag", default=False, action="store_true",
            help="add tag (/1, /2) to the read name")
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(p.print_help())

    tag = opts.tag

    pairsfastqs = args[:-1]
    base = args[-1]
    afastq = base + ".1.fastq"
    bfastq = base + ".2.fastq"
    assert not op.exists(afastq)

    afw = open(afastq, "w")
    bfw = open(bfastq, "w")

    for pairsfastq in pairsfastqs:
        assert op.exists(pairsfastq)
        it = iter_fastq(pairsfastq)
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

    for f in (afw, bfw):
        f.close()

    logging.debug("reads unpaired into `{0}` and `{1}`".\
            format(afastq, bfastq))


def pairinplace(args):
    """
    %prog pairinplace bulk.fastq

    Pair up the records in bulk.fastq by comparing the names for adjancent
    records. If they match, print to bulk.pairs.fastq, else print to
    bulk.frags.fastq.
    """
    from jcvi.utils.iter import pairwise

    p = OptionParser(pairinplace.__doc__)
    p.add_option("-r", dest="rclip", default=1, type="int",
            help="pair ID is derived from rstrip N chars [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    base = op.basename(fastqfile).split(".")[0]

    frags = base + ".frags.fastq"
    pairs = base + ".pairs.fastq"
    if fastqfile.endswith(".gz"):
        frags += ".gz"
        pairs += ".gz"

    fragsfw = must_open(frags, "w")
    pairsfw = must_open(pairs, "w")

    N = opts.rclip
    strip_name = lambda x: x[:-N] if N else str

    fh_iter = iter_fastq(fastqfile, key=strip_name)
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

    logging.debug("Reads paired into `%s` and `%s`" % (pairs, frags))


def pair(args):
    """
    %prog pair 1.fastq 2.fastq

    Pair up the records in 1.fastq and 2.fastq, pairs are indicated by trailing
    "/1" and "/2". If using raw sequences, this is trivial, since we can just
    iterate one at a time for both files; however if two files do not match,
    (e.g. due to trimming), we need a third fastq that provides the order. Two
    output files will be automatically written, one `frags.fastq` and
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
    frags = base + ".frags.fastq"
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
        totalsize = getfilesize(ref)
    else:
        totalsize = getfilesize(afastq)

    from jcvi.apps.console import ProgressBar

    bar = ProgressBar(maxval=totalsize).start()
    strip_name = lambda x: x.rsplit("/", 1)[0]

    ah_iter = iter_fastq(afastq, offset=offset, key=strip_name)
    bh_iter = iter_fastq(bfastq, offset=offset, key=strip_name)

    a = ah_iter.next()
    b = bh_iter.next()

    fragsfw = open(frags, "w")
    pairsfw = open(pairs, "w")

    if ref:
        for r in iter_fastq(ref, offset=0, key=strip_name):
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
