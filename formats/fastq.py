#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Processing fastq files
"""

import os.path as op
import sys
import re
import logging

from itertools import islice

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from jcvi.formats.fasta import must_open, rc
from jcvi.formats.base import DictFile
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


qual_offset = lambda x: 33 if x == "sanger" else 64


class FastqLite (object):
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual

    def __str__(self):
        return "\n".join((self.name, self.seq, "+", self.qual))

    def rc(self):
        self.seq = rc(self.seq)
        self.qual = self.qual[::-1]


class FastqRecord (object):
    def __init__(self, fh, offset=0, key=None):
        self.name = self.header = fh.readline()
        if not self.name:
            return
        self.name = self.name.split()[0]
        self.seq = fh.readline().rstrip()
        self.l3 = fh.readline().rstrip()
        self.qual = fh.readline().rstrip()
        if offset != 0:
            self.qual = "".join(chr(ord(x) + offset) for x in self.qual)
        self.length = len(self.seq)
        assert self.length == len(self.qual), \
                "length mismatch: seq(%s) and qual(%s)" % (self.seq, self.qual)
        if key:
            self.name = key(self.name)

    def __str__(self):
        return "\n".join((self.name, self.seq, "+", self.qual))

    def __len__(self):
        return self.length

    @property
    def quality(self):
        return [ord(x) for x in self.qual]


class FastqHeader(object):

    def __init__(self, row, paired=False):
        header = row.strip().split(" ")
        h = header[0].split(":")
        self.instrument = h[0]
        if len(header) == 2 and header[1].find(":"):
            self.dialect = ">=1.8"  # Illumina Casava 1.8+ format

            self.runId = int(h[1])
            self.flowcellId = h[2]
            self.laneNum = int(h[3])
            self.tileNum = int(h[4])
            self.xPos = int(h[5])
            self.yPos = h[6]
            if re.search("/", self.yPos):
                self.paired = True
                self.yPos, self.readNum = self.yPos.split("/")

            a = header[1].split(":")
            self.readNum = int(a[0])
            self.isFiltered = a[1]
            self.controlNum = int(a[2])
            self.barcode = a[3] if a[3] else 0
        else:
            self.dialect = "<1.8"   # Old Illumina Casava format (< 1.8)
            self.laneNum = int(h[1])
            self.tileNum = int(h[2])
            self.xPos = int(h[3])
            self.yPos = h[4]
            self.paired = False
            m = re.search(r"(\d+)(#\S+)\/(\d+)", self.yPos)
            if m:
                self.paired = True
                self.yPos, self.multiplexId, self.readNum = \
                        m.group(1), m.group(2), m.group(3)

        self.paired = paired

    def __str__(self):
        if self.dialect == ">=1.8":
            yPos = "{0}/{1}".format(self.yPos, self.readNum) if self.paired \
                    else self.yPos

            h0 = ":".join(str(x) for x in (self.instrument, self.runId, \
                    self.flowcellId, self.laneNum, self.tileNum, \
                    self.xPos, yPos))
            h1 = ":".join(str(x) for x in (self.readNum, self.isFiltered, \
                    self.controlNum, self.barcode))

            return "{0} {1}".format(h0, h1)
        else:
            yPos = "{0}#{1}/{2}".format(self.yPos, self.multiplexId, \
                    self.readNum) if self.paired else self.yPos

            return ":".join(str(x) for x in (self.instrument, self.laneNum, \
                    self.tileNum, self.xPos, self.yPos))

    @property
    def illumina_old(self):
        header = ":".join(str(x) for x in (self.instrument, self.laneNum, \
                self.tileNum, self.xPos, self.yPos + "#" + self.barcode))
        header = header + "/" + str(self.readNum)

        return header

    @property
    def is_new_fmt(self):
        return True if self.dialect == ">=1.8" else None


def pairspf(pp):
    return op.basename(op.commonprefix(pp).rstrip("._-"))


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
        ('catread', 'cat pairs together (reverse of splitread)'),
        ('pairinplace', 'collect pairs by checking adjacent ids'),
        ('convert', 'convert between illumina and sanger offset'),
        ('first', 'get first N reads from file'),
        ('filter', 'filter to get high qv reads'),
        ('trim', 'trim reads using fastx_trimmer'),
        ('some', 'select a subset of fastq reads'),
        ('guessoffset', 'guess the quality offset of the fastq records'),
        ('readlen', 'calculate read length'),
        ('format', 'format fastq file, convert header from casava 1.8+ to older format'),
        ('fasta', 'convert fastq to fasta and qual file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def calc_readlen(f, first):
    from jcvi.utils.cbook import SummaryStats

    L = []
    ai = iter_fastq(f)
    rec = ai.next()
    while rec:
        L.append(rec.length)
        if len(L) > first:
            break
        rec = ai.next()
    s = SummaryStats(L)

    return s


def readlen(args):
    """
    %prog readlen fastqfile

    Calculate read length, will only try the first N reads. Output min, max, and
    avg for each file.
    """
    p = OptionParser(readlen.__doc__)
    p.set_firstN()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    f, = args
    s = calc_readlen(f, opts.firstN)
    print "\t".join(str(x) for x in (f, s.min, s.max, s.mean))
    return int(s.max)


def fasta(args):
    """
    %prog fasta fastqfile

    Convert fastq to fasta and qual file.
    """
    p = OptionParser(fasta.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    pf = fastqfile.rsplit(".", 1)[0]
    fastafile, qualfile = pf + ".fasta", pf + ".qual"
    SeqIO.convert(fastqfile, "fastq", fastafile, "fasta")
    SeqIO.convert(fastqfile, "fastq", qualfile, "qual")

    return fastafile, qualfile


def first(args):
    """
    %prog first N fastqfile(s)

    Get first N reads from file.
    """
    from jcvi.apps.base import need_update

    p = OptionParser(first.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    N = int(args[0])
    nlines = N * 4
    fastqfiles = args[1:]
    fastqfile = fastqfiles[0]
    outfile = opts.outfile
    if not need_update(fastqfiles, outfile):
        logging.debug("File `{0}` exists. Will not overwrite.".format(outfile))
        return

    gz = fastqfile.endswith(".gz")
    for fastqfile in fastqfiles:
        if gz:
            cmd = "zcat {0} | head -n {1}".format(fastqfile, nlines)
        else:
            cmd = "head -n {0} {1}".format(nlines, fastqfile)

        sh(cmd, outfile=opts.outfile, append=True)


def FastqPairedIterator(read1, read2):
    if read1 == read2:
        p1fp = p2fp = must_open(read1)
    else:
        p1fp = must_open(read1)
        p2fp = must_open(read2)

    return p1fp, p2fp


def isHighQv(qs, qvchar, pct=90):
    cutoff = len(qs) * pct / 100
    highs = sum(1 for x in qs if x >= qvchar)
    return highs >= cutoff


def filter(args):
    """
    %prog filter paired.fastq

    Filter to get high qv reads. Use interleaved format (one file) or paired
    format (two files) to filter on paired reads.
    """
    p = OptionParser(filter.__doc__)
    p.add_option("-q", dest="qv", default=20, type="int",
                 help="Minimum quality score to keep [default: %default]")
    p.add_option("-p", dest="pct", default=95, type="int",
                 help="Minimum percent of bases that have [-q] quality "\
                 "[default: %default]")

    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    if len(args) == 1:
        r1 = r2 = args[0]
    else:
        r1, r2 = args

    qv = opts.qv
    pct = opts.pct

    offset = guessoffset([r1])
    qvchar = chr(offset + qv)
    logging.debug("Call base qv >= {0} as good.".format(qvchar))
    outfile = r1.rsplit(".", 1)[0] + ".q{0}.paired.fastq".format(qv)
    fw = open(outfile, "w")

    p1fp, p2fp = FastqPairedIterator(r1, r2)
    while True:
        a = list(islice(p1fp, 4))
        if not a:
            break

        b = list(islice(p2fp, 4))
        q1 = a[-1].rstrip()
        q2 = b[-1].rstrip()

        if isHighQv(q1, qvchar, pct=pct) and isHighQv(q2, qvchar, pct=pct):
            fw.writelines(a)
            fw.writelines(b)


def checkShuffleSizes(p1, p2, pairsfastq, extra=0):
    from jcvi.apps.base import getfilesize

    pairssize = getfilesize(pairsfastq)
    p1size = getfilesize(p1)
    p2size = getfilesize(p2)
    assert pairssize == p1size + p2size + extra, \
          "The sizes do not add up: {0} + {1} + {2} != {3}".\
          format(p1size, p2size, extra, pairssize)


def shuffle(args):
    """
    %prog shuffle p1.fastq p2.fastq

    Shuffle pairs into interleaved format.
    """
    p = OptionParser(shuffle.__doc__)
    p.set_tag()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    p1, p2 = args
    pairsfastq = pairspf((p1, p2)) + ".fastq"
    tag = opts.tag

    p1fp = must_open(p1)
    p2fp = must_open(p2)
    pairsfw = must_open(pairsfastq, "w")
    nreads = 0
    while True:
        a = list(islice(p1fp, 4))
        if not a:
            break

        b = list(islice(p2fp, 4))
        if tag:
            name = a[0].rstrip()
            a[0] = name + "/1\n"
            b[0] = name + "/2\n"

        pairsfw.writelines(a)
        pairsfw.writelines(b)
        nreads += 2

    pairsfw.close()
    extra = nreads * 2 if tag else 0
    checkShuffleSizes(p1, p2, pairsfastq, extra=extra)

    logging.debug("File `{0}` verified after writing {1} reads.".\
                     format(pairsfastq, nreads))
    return pairsfastq


def split(args):
    """
    %prog split pairs.fastq

    Split shuffled pairs into `.1.fastq` and `.2.fastq`, using `sed`. Can work
    on gzipped file.

    <http://seqanswers.com/forums/showthread.php?t=13776>
    """
    from jcvi.apps.grid import Jobs

    p = OptionParser(split.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pairsfastq, = args
    gz = pairsfastq.endswith(".gz")
    pf = pairsfastq.replace(".gz", "").rsplit(".", 1)[0]
    p1 = pf + ".1.fastq"
    p2 = pf + ".2.fastq"

    cmd = "zcat" if gz else "cat"
    p1cmd = cmd + " {0} | sed -ne '1~8{{N;N;N;p}}'".format(pairsfastq)
    p2cmd = cmd + " {0} | sed -ne '5~8{{N;N;N;p}}'".format(pairsfastq)

    if gz:
        p1cmd += " | gzip"
        p2cmd += " | gzip"
        p1 += ".gz"
        p2 += ".gz"

    p1cmd += " > " + p1
    p2cmd += " > " + p2

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


def format(args):
    """
    %prog format fastqfile

    Format FASTQ file. Currently provides option to convert FASTQ header from
    Illumina Casava 1.8+ format to the older format
    """
    p = OptionParser(format.__doc__)

    p.add_option("--old_header", default=False, action="store_true",
                help="Convert header format from illumina new (1.8+) to older format" +
                " [default: %default]")
    p.set_tag()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    ai = iter_fastq(fastqfile)
    rec = ai.next()
    while rec:
        h = FastqHeader(rec.header, paired=opts.tag)
        if opts.old_header:
            if h.is_new_fmt:
                rec.name = h.illumina_old
            else:
                logging.error("Error: Input fastq header not in Illumina Casava" +
                            " 1.8+ format")
                sys.exit()
        else:
            rec.name = str(h)
        print rec
        rec = ai.next()


def some(args):
    """
    %prog some idsfile afastq [bfastq]

    Select a subset of the reads with ids present in the idsfile.
    `bfastq` is optional (only if reads are paired)
    """
    p = OptionParser(some.__doc__)
    opts, args = p.parse_args(args)

    if len(args) not in (2, 3):
        sys.exit(not p.print_help())

    idsfile, afastq, = args[:2]
    bfastq = args[2] if len(args) == 3 else None

    ids = DictFile(idsfile, valuepos=None)

    ai = iter_fastq(open(afastq))
    arec = ai.next()
    if bfastq:
        bi = iter_fastq(open(bfastq))
        brec = bi.next()

    while arec:
        if arec.name[1:] in ids:
            print arec
            if bfastq:
                print brec

        arec = ai.next()
        if bfastq:
            brec = bi.next()


def trim(args):
    """
    %prog trim fastqfile

    Wraps `fastx_trimmer` to trim from begin or end of reads.
    """
    p = OptionParser(trim.__doc__)
    p.add_option("-f", dest="first", default=0, type="int",
            help="First base to keep. Default is 1.")
    p.add_option("-l", dest="last", default=0, type="int",
            help="Last base to keep. Default is entire read.")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    obfastqfile = op.basename(fastqfile)
    fq = obfastqfile.rsplit(".", 1)[0] + ".ntrimmed.fastq"
    if fastqfile.endswith(".gz"):
        fq = obfastqfile.rsplit(".", 2)[0] + ".ntrimmed.fastq.gz"

    cmd = "fastx_trimmer -Q33 "
    if opts.first:
        cmd += "-f {0.first} ".format(opts)
    if opts.last:
        cmd += "-l {0.last} ".format(opts)

    sh(cmd, infile=fastqfile, outfile=fq)


def catread(args):
    """
    %prog catread fastqfile1 fastqfile2

    Concatenate paired end reads into one. Useful for example to do single-end
    mapping and perform filtering on the whole read pair level.
    """
    p = OptionParser(catread.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    r1, r2 = args
    p1fp, p2fp = FastqPairedIterator(r1, r2)
    outfile = pairspf((r1, r2)) + ".cat.fastq"
    fw = must_open(outfile, "w")
    while True:
        a = list(islice(p1fp, 4))
        if not a:
            break
        atitle, aseq, _, aqual = a
        btitle, bseq, _, bqual = list(islice(p2fp, 4))
        print >> fw, "\n".join((atitle.strip(), aseq.strip() + bseq.strip(), \
                                "+", aqual.strip() + bqual.strip()))


def splitread(args):
    """
    %prog splitread fastqfile

    Split fastqfile into two read fastqfiles, cut in the middle.
    """
    p = OptionParser(splitread.__doc__)
    p.add_option("-n", dest="n", default=76, type="int",
            help="Split at N-th base position [default: %default]")
    p.add_option("--rc", default=False, action="store_true",
            help="Reverse complement second read [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pairsfastq, = args

    base = op.basename(pairsfastq).split(".")[0]
    fq1 = base + ".1.fastq"
    fq2 = base + ".2.fastq"
    fw1 = must_open(fq1, "w")
    fw2 = must_open(fq2, "w")

    fp = must_open(pairsfastq)
    n = opts.n

    for name, seq, qual in FastqGeneralIterator(fp):

        name = "@" + name
        rec1 = FastqLite(name, seq[:n], qual[:n])
        rec2 = FastqLite(name, seq[n:], qual[n:])
        if opts.rc:
            rec2.rc()

        print >> fw1, rec1
        print >> fw2, rec2

    logging.debug("Reads split into `{0},{1}`".format(fq1, fq2))
    fw1.close()
    fw2.close()


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
    p.set_phred()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    infastq, outfastq = args
    phred = opts.phred or str(guessoffset([infastq]))
    fin = "illumina" if phred == "64" else "sanger"
    fout = "sanger" if phred == "64" else "illumina"

    seqret = "seqret"
    if infastq.endswith(".gz"):
        cmd = "zcat {0} | ".format(infastq)
        cmd += seqret + " fastq-{0}::stdin fastq-{1}::stdout".\
                format(fin, fout)
    else:
        cmd = seqret + " fastq-{0}::{1} fastq-{2}::stdout".\
                format(fin, infastq, fout)

    sh(cmd, outfile=outfastq)

    return outfastq


def pairinplace(args):
    """
    %prog pairinplace bulk.fastq

    Pair up the records in bulk.fastq by comparing the names for adjancent
    records. If they match, print to bulk.pairs.fastq, else print to
    bulk.frags.fastq.
    """
    from jcvi.utils.iter import pairwise

    p = OptionParser(pairinplace.__doc__)
    p.set_rclip()
    p.set_tag()
    p.add_option("--base",
                help="Base name for the output files [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    base = opts.base or op.basename(fastqfile).split(".")[0]

    frags = base + ".frags.fastq"
    pairs = base + ".pairs.fastq"
    if fastqfile.endswith(".gz"):
        frags += ".gz"
        pairs += ".gz"

    fragsfw = must_open(frags, "w")
    pairsfw = must_open(pairs, "w")

    N = opts.rclip
    tag = opts.tag
    strip_name = (lambda x: x[:-N]) if N else None

    fh_iter = iter_fastq(fastqfile, key=strip_name)
    skipflag = False  # controls the iterator skip
    for a, b in pairwise(fh_iter):
        if b is None:  # hit the eof
            break

        if skipflag:
            skipflag = False
            continue

        if a.name == b.name:
            if tag:
                a.name += "/1"
                b.name += "/2"
            print >> pairsfw, a
            print >> pairsfw, b
            skipflag = True
        else:
            print >> fragsfw, a

    # don't forget the last one, when b is None
    if not skipflag:
        print >> fragsfw, a

    logging.debug("Reads paired into `%s` and `%s`" % (pairs, frags))
    return pairs


if __name__ == '__main__':
    main()
