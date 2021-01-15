"""
Wrapper for biopython Fasta, add option to parse sequence headers
"""
from __future__ import print_function

import re
import sys
import os
import os.path as op
import shutil
import logging
import string
import hashlib

from itertools import groupby, zip_longest
from more_itertools import grouper, pairwise

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

from jcvi.formats.base import BaseFile, DictFile, must_open
from jcvi.formats.bed import Bed
from jcvi.utils.cbook import percentage
from jcvi.utils.console import printf
from jcvi.utils.table import write_csv
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update


class Fasta(BaseFile, dict):
    def __init__(self, filename, index=False, key_function=None, lazy=False):
        super(Fasta, self).__init__(filename)
        self.key_function = key_function

        if lazy:  # do not incur the overhead
            return

        if index:
            self.index = SeqIO.index(filename, "fasta", key_function=key_function)
        else:
            # SeqIO.to_dict expects a different key_function that operates on
            # the SeqRecord instead of the raw string
            _key_function = (
                (lambda rec: key_function(rec.description)) if key_function else None
            )
            self.index = SeqIO.to_dict(
                SeqIO.parse(must_open(filename), "fasta"), key_function=_key_function
            )

    def _key_function(self, key):
        return self.key_function(key) if self.key_function else key

    def __len__(self):
        return len(self.index)

    def __contains__(self, key):
        key = self._key_function(key)
        return key in self.index

    def __getitem__(self, key):
        key = self._key_function(key)
        rec = self.index[key]
        return rec

    def keys(self):
        return self.index.keys()

    def iterkeys(self):
        for k in self.index.keys():
            yield k

    def iteritems(self):
        for k in self.iterkeys():
            yield k, self[k]

    def itersizes(self):
        for k in self.iterkeys():
            yield k, len(self[k])

    def iteritems_ordered(self):
        for rec in SeqIO.parse(must_open(self.filename), "fasta"):
            yield rec.name, rec

    def iterdescriptions_ordered(self):
        for k, rec in self.iteritems_ordered():
            yield rec.description, rec

    def iterkeys_ordered(self):
        for k, rec in self.iteritems_ordered():
            yield k

    def itersizes_ordered(self):
        for k, rec in self.iteritems_ordered():
            yield k, len(rec)

    def tostring(self):
        d = {}
        for k, rec in self.iteritems():
            d[k] = str(rec.seq)
        return d

    @property
    def totalsize(self):
        return sum(size for k, size in self.itersizes())

    @classmethod
    def subseq(cls, fasta, start=None, stop=None, strand=None):
        """
        Take Bio.SeqRecord and slice "start:stop" from it, does proper
        index and error handling
        """
        start = start - 1 if start is not None else 0
        stop = stop if stop is not None else len(fasta)

        if start < 0:
            msg = "start ({0}) must > 0 of `{1}`. Reset to 1".format(
                start + 1, fasta.id
            )
            logging.error(msg)
            start = 0

        if stop > len(fasta):
            msg = "stop ({0}) must be <= length of `{1}` ({2}). Reset to {2}.".format(
                stop, fasta.id, len(fasta)
            )
            logging.error(msg)
            stop = len(fasta)

        seq = fasta.seq[start:stop]

        if strand in (-1, "-1", "-"):
            seq = seq.reverse_complement()

        return seq

    def sequence(self, f, asstring=True):
        """
        Emulate brentp's pyfasta/fasta.py sequence() methods

        take a feature and use the start/stop or exon_keys to return
        the sequence from the assocatied fasta file:

        f: a feature
        asstring: if true, return the sequence as a string
                : if false, return as a biopython Seq

        >>> f = Fasta('tests/data/three_chrs.fasta')
        >>> f.sequence({'start':1, 'stop':2, 'strand':1, 'chr': 'chr1'})
        'AC'
        >>> f.sequence({'start':1, 'stop':2, 'strand': -1, 'chr': 'chr1'})
        'GT'
        """

        assert "chr" in f, "`chr` field required"
        name = f["chr"]

        assert name in self, "feature: %s not in `%s`" % (f, self.filename)

        fasta = self[f["chr"]]

        seq = Fasta.subseq(fasta, f.get("start"), f.get("stop"), f.get("strand"))

        if asstring:
            return str(seq)

        return seq


"""
Class derived from https://gist.github.com/933737
Original code written by David Winter (https://github.com/dwinter)

Code writted to answer this challenge at Biostar:
http://biostar.stackexchange.com/questions/5902/

(Code includes improvements from Brad Chapman)
"""


class ORFFinder(object):
    """Find the longest ORF in a given sequence
    "seq" is a string, if "start" is not provided any codon can be the start of
    and ORF. If muliple ORFs have the longest length the first one encountered
    is printed
    """

    def __init__(self, seq, start=[], stop=["TAG", "TAA", "TGA"]):
        self.seq = str(seq).upper()
        self.start = start
        self.stop = stop
        # strand, frame, start, end, length; coordinates are 1-based
        self.result = ["+", 0, 0, 0, 0]
        self.longest = 0
        self.size = len(seq)

    def __str__(self):
        # Format similar to getorf
        strand, frame, start, end, length = self.result
        start += 1  # 1-based coordinates
        if strand == "-":
            start, end = end, start
        return "[{0} - {1}]".format(start, end)

    @property
    def info(self):
        strand, frame, start, end, length = self.result
        return "\t".join(str(x) for x in (strand, frame, start, end))

    def codons(self, frame):
        """A generator that yields DNA in one codon blocks
        "frame" counts for 0. This function yields a tuple (triplet, index) with
        index relative to the original DNA sequence
        """
        start = frame
        while start + 3 <= self.size:
            yield self.sequence[start : start + 3], start
            start += 3

    def scan_sequence(self, frame, direction):
        """ Search in one reading frame """
        orf_start = None
        for c, index in self.codons(frame):
            if (
                c not in self.stop
                and (c in self.start or not self.start)
                and orf_start is None
            ):
                orf_start = index
            elif c in self.stop and orf_start is not None:
                self._update_longest(orf_start, index + 3, direction, frame)
                orf_start = None

        if orf_start is not None:
            self._update_longest(orf_start, index + 3, direction, frame)

    def _update_longest(self, orf_start, index, direction, frame):
        orf_end = index
        L = orf_end - orf_start
        if L > self.longest:
            self.longest = L
            self.result = [direction, frame, orf_start, orf_end, L]

    def get_longest_orf(self):
        dirs = ("+", "-")
        for direction in dirs:
            self.sequence = self.seq
            if direction == "-":
                self.sequence = rc(self.sequence)
            for frame in range(3):
                self.scan_sequence(frame, direction)

        strand, frame, start, end, length = self.result
        size = self.size
        if strand == "-":
            start, end = size - end, size - start
            self.result[2:4] = start, end

        assert start <= end, self.result
        if start == end:
            return "N"

        orf = self.seq[start:end]
        if strand == "-":
            orf = rc(orf)

        assert len(orf) % 3 == 0

        return orf


class SequenceInfo(object):
    """
    Emulate output from `sequence_info`:

    File                           SUBAC32.contigs.fasta

    Number of sequences                    80

    Residue counts:
      Number of A's                     66266   31.36 %
      Number of C's                     40032   18.95 %
      Number of G's                     39145   18.53 %
      Number of T's                     65799   31.14 %
      Number of N's                        58    0.03 %
      Total                            211300

    Sequence lengths:
      Minimum                             242
      Maximum                            8398
      Average                            2641.25
      N50                                4791
    """

    def __init__(self, filename, gapstats=False):
        from jcvi.utils.cbook import SummaryStats
        from jcvi.assembly.base import calculate_A50

        f = Fasta(filename)
        self.filename = filename
        self.header = "File|#_seqs|#_reals|#_Ns|Total|Min|Max|N50".split("|")
        if gapstats:
            self.header += ["Gaps"]
        self.nseqs = len(f)
        sizes = []
        gaps = []
        na = nc = ng = nt = 0
        for k, s in f.iteritems():
            s = str(s.seq).upper()
            sizes.append(len(s))
            na += s.count("A")
            nc += s.count("C")
            ng += s.count("G")
            nt += s.count("T")
            if gapstats:
                gaps += list(self.iter_gap_len(s))
        self.real = real = na + nc + ng + nt
        s = SummaryStats(sizes)
        self.sum = s.sum
        if gapstats:
            self.gaps = len(gaps)
        self.nn = self.sum - real
        a50, l50, nn50 = calculate_A50(sizes)
        self.min = s.min
        self.max = s.max
        self.mean = int(s.mean)
        self.n50 = l50
        self.data = [
            self.filename,
            self.nseqs,
            self.real,
            self.nn,
            self.sum,
            self.min,
            self.max,
            self.n50,
        ]
        if gapstats:
            self.data += [self.gaps]
        assert len(self.header) == len(self.data)

    def iter_gap_len(self, seq, mingap=10):
        for gap, seq in groupby(seq, lambda x: x == "N"):
            if not gap:
                continue
            gap_len = len(list(seq))
            if gap_len >= mingap:
                yield len(list(seq))


def rc(s):
    _complement = str.maketrans("ATCGatcgNnXx", "TAGCtagcNnXx")
    cs = s.translate(_complement)
    return cs[::-1]


def main():

    actions = (
        (
            "extract",
            "given fasta file and seq id, retrieve the sequence " + "in fasta format",
        ),
        ("longestorf", "find longest orf for CDS fasta"),
        ("translate", "translate CDS to proteins"),
        ("info", "run `sequence_info` on fasta files"),
        ("summary", "report the real no of bases and N's in fasta files"),
        ("uniq", "remove records that are the same"),
        ("ids", "generate a list of headers"),
        (
            "format",
            "trim accession id to the first space or switch id "
            + "based on 2-column mapping file",
        ),
        ("pool", "pool a bunch of fastafiles together and add prefix"),
        ("random", "randomly take some records"),
        ("simulate", "simulate random fasta file for testing"),
        ("diff", "check if two fasta records contain same information"),
        ("identical", "given 2 fasta files, find all exactly identical records"),
        ("trim", "given a cross_match screened fasta, trim the sequence"),
        ("trimsplit", "split sequences at lower-cased letters"),
        ("sort", "sort the records by IDs, sizes, etc."),
        ("filter", "filter the records by size"),
        ("pair", "sort paired reads to .pairs, rest to .fragments"),
        (
            "pairinplace",
            "starting from fragment.fasta, find if "
            + "adjacent records can form pairs",
        ),
        ("fastq", "combine fasta and qual to create fastq file"),
        ("tidy", "normalize gap sizes and remove small components in fasta"),
        ("sequin", "generate a gapped fasta file for sequin submission"),
        ("gaps", "print out a list of gap sizes within sequences"),
        ("join", "concatenate a list of seqs and add gaps in between"),
        (
            "some",
            "include or exclude a list of records (also performs on "
            + ".qual file if available)",
        ),
        ("qual", "generate dummy .qual file based on FASTA file"),
        ("clean", "remove irregular chars in FASTA seqs"),
        ("ispcr", "reformat paired primers into isPcr query format"),
        ("fromtab", "convert 2-column sequence file to FASTA format"),
        ("gc", "plot G+C content distribution"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def simulate_one(fw, name, size):
    """
    Simulate a random sequence with name and size
    """
    from random import choice

    seq = Seq("".join(choice("ACGT") for _ in range(size)))
    s = SeqRecord(seq, id=name, description="Fake sequence")
    SeqIO.write([s], fw, "fasta")


def simulate(args):
    """
    %prog simulate idsfile

    Simulate random FASTA file based on idsfile, which is a two-column
    tab-separated file with sequence name and size.
    """
    p = OptionParser(simulate.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (idsfile,) = args
    fp = open(idsfile)
    fw = must_open(opts.outfile, "w")
    for row in fp:
        name, size = row.split()
        size = int(size)
        simulate_one(fw, name, size)
    fp.close()


def gc(args):
    """
    %prog gc fastafile

    Plot G+C content distribution.
    """
    p = OptionParser(gc.__doc__)
    p.add_option("--binsize", default=500, type="int", help="Bin size to use")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    binsize = opts.binsize
    allbins = []
    for name, seq in parse_fasta(fastafile):
        for i in range(len(seq) / binsize):
            atcnt = gccnt = 0
            for c in seq[i * binsize : (i + 1) * binsize].upper():
                if c in "AT":
                    atcnt += 1
                elif c in "GC":
                    gccnt += 1
            totalcnt = atcnt + gccnt
            if totalcnt == 0:
                continue
            gcpct = gccnt * 100 / totalcnt
            allbins.append(gcpct)

    from jcvi.graphics.base import asciiplot
    from collections import Counter

    title = "Total number of bins={}".format(len(allbins))
    c = Counter(allbins)
    x, y = zip(*sorted(c.items()))
    asciiplot(x, y, title=title)


def trimsplit(args):
    """
    %prog trimsplit fastafile

    Split sequences at lower-cased letters and stretch of Ns. This is useful
    at cleaning up the low quality bases for the QUIVER output.
    """
    from jcvi.utils.cbook import SummaryStats

    p = OptionParser(trimsplit.__doc__)
    p.add_option(
        "--minlength", default=1000, type="int", help="Min length of contigs to keep"
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    minlength = opts.minlength

    fw = must_open(fastafile.rsplit(".", 1)[0] + ".split.fasta", "w")
    ntotal = 0
    removed = []
    Ns = []
    for name, seq in parse_fasta(fastafile):
        stretches = []
        ntotal += len(seq)
        for lower, stretch in groupby(seq, key=lambda x: x.islower()):
            stretch = "".join(stretch)
            if lower or len(stretch) < minlength:
                removed.append(len(stretch))
                continue
            for isN, s in groupby(stretch, key=lambda x: x in "Nn"):
                s = "".join(s)
                if isN or len(s) < minlength:
                    Ns.append(len(s))
                    continue
                stretches.append(s)
        for i, seq in enumerate(stretches):
            id = "{0}_{1}".format(name.split("|")[0], i)
            s = SeqRecord(Seq(seq), id=id, description="")
            SeqIO.write([s], fw, "fasta")
    fw.close()

    # Reporting
    if removed:
        logging.debug(
            "Total bases removed: {0}".format(percentage(sum(removed), ntotal))
        )
        print(SummaryStats(removed), file=sys.stderr)
    if Ns:
        logging.debug("Total Ns removed: {0}".format(percentage(sum(Ns), ntotal)))
        print(SummaryStats(Ns), file=sys.stderr)


def qual(args):
    """
    %prog qual fastafile

    Generate dummy .qual file based on FASTA file.
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(qual.__doc__)
    p.add_option(
        "--qv", default=31, type="int", help="Dummy qv score for extended bases"
    )
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    sizes = Sizes(fastafile)
    qvchar = str(opts.qv)
    fw = must_open(opts.outfile, "w")
    total = 0
    for s, slen in sizes.iter_sizes():
        print(">" + s, file=fw)
        print(" ".join([qvchar] * slen), file=fw)
        total += 1
    fw.close()
    logging.debug("Written {0} records in `{1}`.".format(total, opts.outfile))


def info(args):
    """
    %prog info *.fasta

    Run `sequence_info` on FASTA files. Generate a report per file.
    """
    p = OptionParser(info.__doc__)
    p.add_option(
        "--gaps", default=False, action="store_true", help="Count number of gaps"
    )
    p.set_table()
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    fastafiles = args
    data = []
    for f in fastafiles:
        s = SequenceInfo(f, gapstats=opts.gaps)
        data.append(s.data)
    write_csv(s.header, data, sep=opts.sep, filename=opts.outfile, align=opts.align)


def fromtab(args):
    """
    %prog fromtab tabfile fastafile

    Convert 2-column sequence file to FASTA format. One usage for this is to
    generatea `adapters.fasta` for TRIMMOMATIC.
    """
    p = OptionParser(fromtab.__doc__)
    p.set_sep(sep=None)
    p.add_option(
        "--noheader", default=False, action="store_true", help="Ignore first line"
    )
    p.add_option("--replace", help="Replace spaces in name to char")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    tabfile, fastafile = args
    sep = opts.sep
    replace = opts.replace
    fp = must_open(tabfile)
    fw = must_open(fastafile, "w")
    nseq = 0
    if opts.noheader:
        next(fp)
    for row in fp:
        row = row.strip()
        if not row or row[0] == "#":
            continue

        name, seq = row.rsplit(sep, 1)
        if replace:
            name = name.replace(" ", replace)
        print(">{0}\n{1}".format(name, seq), file=fw)
        nseq += 1
    fw.close()

    logging.debug("A total of {0} sequences written to `{1}`.".format(nseq, fastafile))


def longestorf(args):
    """
    %prog longestorf fastafile

    Find longest ORF for each sequence in fastafile.
    """
    p = OptionParser(longestorf.__doc__)
    p.add_option("--ids", action="store_true", help="Generate table with ORF info")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    pf = fastafile.rsplit(".", 1)[0]
    orffile = pf + ".orf.fasta"
    idsfile = None
    if opts.ids:
        idsfile = pf + ".orf.ids"
        fwids = open(idsfile, "w")

    f = Fasta(fastafile, lazy=True)
    fw = must_open(orffile, "w")
    before, after = 0, 0
    for name, rec in f.iteritems_ordered():
        cds = rec.seq
        before += len(cds)
        # Try all six frames
        orf = ORFFinder(cds)
        lorf = orf.get_longest_orf()
        newcds = Seq(lorf)
        after += len(newcds)
        newrec = SeqRecord(newcds, id=name, description=rec.description)
        SeqIO.write([newrec], fw, "fasta")
        if idsfile:
            print("\t".join((name, orf.info)), file=fwids)

    fw.close()
    if idsfile:
        fwids.close()

    logging.debug(
        "Longest ORFs written to `{0}` ({1}).".format(
            orffile, percentage(after, before)
        )
    )

    return orffile


def ispcr(args):
    """
    %prog ispcr fastafile

    Reformat paired primers into isPcr query format, which is three column
    format: name, forward, reverse
    """
    p = OptionParser(ispcr.__doc__)
    p.add_option(
        "-r",
        dest="rclip",
        default=1,
        type="int",
        help="pair ID is derived from rstrip N chars",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    ispcrfile = fastafile + ".isPcr"
    fw = open(ispcrfile, "w")

    N = opts.rclip
    strip_name = lambda x: x[:-N] if N else str

    npairs = 0
    fastaiter = SeqIO.parse(fastafile, "fasta")
    for a, b in grouper(fastaiter, 2):

        aid, bid = [strip_name(x) for x in (a.id, b.id)]
        assert aid == bid, "Name mismatch {0}".format((aid, bid))

        print("\t".join((aid, str(a.seq), str(b.seq))), file=fw)
        npairs += 1

    fw.close()
    logging.debug("A total of {0} pairs written to `{1}`.".format(npairs, ispcrfile))


def parse_fasta(infile, upper=False):
    """
    parse a fasta-formatted file and returns header
    can be a fasta file that contains multiple records.
    """
    try:
        fp = must_open(infile)
    except:
        fp = infile
    # keep header
    fa_iter = (x[1] for x in groupby(fp, lambda row: row[0] == ">"))
    for header in fa_iter:
        header = next(header)
        if header[0] != ">":
            continue
        # drop '>'
        header = header.strip()[1:]
        # stitch the sequence lines together and make into upper case
        seq = "".join(s.strip() for s in next(fa_iter))
        if upper:
            seq = seq.upper()
        yield header, seq


def iter_clean_fasta(fastafile):
    for header, seq in parse_fasta(fastafile):
        seq = "".join(x for x in seq if x in string.letters or x == "*")
        yield header, seq


def iter_canonical_fasta(fastafile):
    canonical = "ACGTN"
    totalbad = 0
    for header, seq in parse_fasta(fastafile):
        badcounts = sum(1 for x in seq if x not in canonical)
        seq = "".join((x if x in canonical else "N") for x in seq)
        totalbad += badcounts
        yield header, seq

    logging.debug("Total bad char: {0}".format(totalbad))


def fancyprint(fw, seq, width=60, chunk=10):
    assert width % chunk == 0
    nchunks = width / chunk
    seqlen = len(seq)
    maxchar = len(str(seqlen))

    s = ["".join(x) for x in grouper(seq, chunk, fillvalue="")]
    s = [" ".join(x) for x in grouper(s, nchunks, fillvalue="")]
    for a, b in zip(range(1, len(seq), width), s):
        b = b.rstrip()
        a = str(a).rjust(maxchar, " ")
        print("  ".join((a, b)), file=fw)


def clean(args):
    """
    %prog clean fastafile

    Remove irregular chars in FASTA seqs.
    """
    p = OptionParser(clean.__doc__)
    p.add_option(
        "--fancy", default=False, action="store_true", help="Pretty print the sequence"
    )
    p.add_option(
        "--canonical", default=False, action="store_true", help="Use only acgtnACGTN"
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    fw = must_open(opts.outfile, "w")
    if opts.fancy:
        for header, seq in iter_clean_fasta(fastafile):
            print(">" + header, file=fw)
            fancyprint(fw, seq)

        return 0

    iterator = iter_canonical_fasta if opts.canonical else iter_clean_fasta

    for header, seq in iterator(fastafile):
        seq = Seq(seq)
        s = SeqRecord(seq, id=header, description="")
        SeqIO.write([s], fw, "fasta")


def translate(args):
    """
    %prog translate cdsfasta

    Translate CDS to proteins. The tricky thing is that sometimes the CDS
    represents a partial gene, therefore disrupting the frame of the protein.
    Check all three frames to get a valid translation.
    """
    from jcvi.utils.cbook import gene_name

    transl_tables = [str(x) for x in range(1, 25)]
    p = OptionParser(translate.__doc__)
    p.add_option(
        "--ids",
        default=False,
        action="store_true",
        help="Create .ids file with the complete/partial/gaps label",
    )
    p.add_option(
        "--longest",
        default=False,
        action="store_true",
        help="Find the longest ORF from each input CDS",
    )
    p.add_option(
        "--table",
        default=1,
        choices=transl_tables,
        help="Specify translation table to use",
    )
    p.add_option(
        "--strip_names",
        default=False,
        action="store_true",
        help="Strip alternative splicing (e.g. At5g06540.1 -> At5g06540)",
    )
    p.add_option(
        "--unique",
        default=False,
        action="store_true",
        help="Ensure the output FASTA contains unique identifiers",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)
    strip_names = opts.strip_names
    unique = opts.unique

    if len(args) != 1:
        sys.exit(not p.print_help())

    (cdsfasta,) = args
    if opts.longest:
        cdsfasta = longestorf([cdsfasta])

    f = Fasta(cdsfasta, lazy=True)
    outfile = opts.outfile
    fw = must_open(outfile, "w")

    if opts.ids:
        idsfile = cdsfasta.rsplit(".", 1)[0] + ".ids"
        ids = open(idsfile, "w")
    else:
        ids = None

    five_prime_missing = three_prime_missing = 0
    contain_ns = complete = cannot_translate = total = 0

    seen = set()
    grand_total = 0
    for name, rec in f.iteritems_ordered():
        grand_total += 1

        if strip_names:
            name = gene_name(name)

        if unique and name in seen:
            continue

        cds = rec.seq
        cdslen = len(cds)
        peplen = cdslen // 3
        total += 1

        # Try all three frames
        pep = ""
        for i in range(3):
            newcds = cds[i : i + peplen * 3]
            newpep = newcds.translate(table=opts.table)
            if len(newpep.split("*")[0]) > len(pep.split("*")[0]):
                pep = newpep

        labels = []
        if "*" in pep.rstrip("*"):
            logging.error("{0} cannot translate".format(name))
            cannot_translate += 1
            labels.append("cannot_translate")

        contains_start = pep.startswith("M")
        contains_stop = pep.endswith("*")
        contains_ns = "X" in pep
        start_ns = pep.startswith("X")
        end_ns = pep.endswith("X")

        if not contains_start:
            five_prime_missing += 1
            labels.append("five_prime_missing")
        if not contains_stop:
            three_prime_missing += 1
            labels.append("three_prime_missing")
        if contains_ns:
            contain_ns += 1
            labels.append("contain_ns")
        if contains_start and contains_stop:
            complete += 1
            labels.append("complete")
        if start_ns:
            labels.append("start_ns")
        if end_ns:
            labels.append("end_ns")

        if ids:
            print("\t".join((name, ",".join(labels))), file=ids)

        peprec = SeqRecord(pep, id=name, description=rec.description)
        SeqIO.write([peprec], fw, "fasta")
        fw.flush()
        seen.add(name)

    print(
        "Complete gene models: {0}".format(percentage(complete, total)), file=sys.stderr
    )
    print(
        "Missing 5`-end: {0}".format(percentage(five_prime_missing, total)),
        file=sys.stderr,
    )
    print(
        "Missing 3`-end: {0}".format(percentage(three_prime_missing, total)),
        file=sys.stderr,
    )
    print("Contain Ns: {0}".format(percentage(contain_ns, total)), file=sys.stderr)

    if cannot_translate:
        print(
            "Cannot translate: {0}".format(percentage(cannot_translate, total)),
            file=sys.stderr,
        )

    fw.close()

    logging.debug(
        "Total records: {}, Unique records (strip_names={}): {}".format(
            grand_total, strip_names, len(seen)
        )
    )

    return cdsfasta, outfile


def filter(args):
    """
    %prog filter fastafile 100

    Filter the FASTA file to contain records with size >= or <= certain cutoff.
    """
    p = OptionParser(filter.__doc__)
    p.add_option(
        "--less",
        default=False,
        action="store_true",
        help="filter the sizes < certain cutoff [default: >=]",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, cutoff = args
    try:
        cutoff = int(cutoff)
    except ValueError:
        sys.exit(not p.print_help())

    f = Fasta(fastafile, lazy=True)

    fw = must_open(opts.outfile, "w")
    for name, rec in f.iteritems_ordered():

        if opts.less and len(rec) >= cutoff:
            continue

        if (not opts.less) and len(rec) < cutoff:
            continue

        SeqIO.write([rec], fw, "fasta")
        fw.flush()

    return fw.name


def pool(args):
    """
    %prog pool fastafiles > pool.fasta

    Pool a bunch of FASTA files, and add prefix to each record based on
    filenames. File names are simplified to longest unique prefix to avoid
    collisions after getting shortened.
    """
    from jcvi.formats.base import longest_unique_prefix

    p = OptionParser(pool.__doc__)
    p.add_option("--sep", default=".", help="Separator between prefix and name")
    p.add_option(
        "--sequential", default=False, action="store_true", help="Add sequential IDs"
    )
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    for fastafile in args:
        pf = longest_unique_prefix(fastafile, args)
        print(fastafile, "=>", pf, file=sys.stderr)
        prefixopt = "--prefix={0}{1}".format(pf, opts.sep)
        format_args = [fastafile, "stdout", prefixopt]
        if opts.sequential:
            format_args += ["--sequential=replace"]
        format(format_args)


def ids(args):
    """
    %prog ids fastafiles

    Generate the FASTA headers without the '>'.
    """
    p = OptionParser(ids.__doc__)
    p.add_option(
        "--until", default=None, help="Truncate the name and description at words"
    )
    p.add_option(
        "--description",
        default=False,
        action="store_true",
        help="Generate a second column with description",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    until = opts.until
    fw = must_open(opts.outfile, "w")
    for row in must_open(args):
        if row[0] == ">":
            row = row[1:].rstrip()
            if until:
                row = row.split(until)[0]

            atoms = row.split(None, 1)
            if opts.description:
                outrow = "\t".join(atoms)
            else:
                outrow = atoms[0]
            print(outrow, file=fw)

    fw.close()


def sort(args):
    """
    %prog sort fastafile

    Sort a list of sequences and output with sorted IDs, etc.
    """
    p = OptionParser(sort.__doc__)
    p.add_option(
        "--sizes", default=False, action="store_true", help="Sort by decreasing size"
    )

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (fastafile,) = args
    sortedfastafile = fastafile.rsplit(".", 1)[0] + ".sorted.fasta"

    f = Fasta(fastafile, index=False)
    fw = must_open(sortedfastafile, "w")
    if opts.sizes:
        # Sort by decreasing size
        sortlist = sorted(f.itersizes(), key=lambda x: (-x[1], x[0]))
        logging.debug(
            "Sort by size: max: {0}, min: {1}".format(sortlist[0], sortlist[-1])
        )
        sortlist = [x for x, s in sortlist]
    else:
        sortlist = sorted(f.iterkeys())

    for key in sortlist:
        rec = f[key]
        SeqIO.write([rec], fw, "fasta")

    logging.debug("Sorted file written to `{0}`.".format(sortedfastafile))
    fw.close()

    return sortedfastafile


def join(args):
    """
    %prog join fastafile [phasefile]

    Make AGP file for a bunch of sequences, and add gaps between, and then build
    the joined fastafile. This is useful by itself, but with --oo option this
    can convert the .oo (BAMBUS output) into AGP and a joined fasta.

    Phasefile is optional, but must contain two columns - BAC and phase (0, 1, 2, 3).
    """
    from jcvi.formats.agp import OO, Phases, build
    from jcvi.formats.sizes import Sizes

    p = OptionParser(join.__doc__)
    p.add_option("--newid", default=None, help="New sequence ID")
    p.add_option(
        "--gapsize",
        default=100,
        type="int",
        help="Number of N's in between the sequences",
    )
    p.add_option("--gaptype", default="contig", help="Gap type to use in the AGP file")
    p.add_option(
        "--evidence", default="", help="Linkage evidence to report in the AGP file"
    )
    p.add_option("--oo", help="Use .oo file generated by bambus")
    opts, args = p.parse_args(args)

    nargs = len(args)
    if nargs not in (1, 2):
        sys.exit(not p.print_help())

    if nargs == 2:
        fastafile, phasefile = args
        phases = DictFile(phasefile)
        phases = dict((a, Phases[int(b)]) for a, b in phases.items())
    else:
        (fastafile,) = args
        phases = {}

    sizes = Sizes(fastafile)
    prefix = fastafile.rsplit(".", 1)[0]
    agpfile = prefix + ".agp"
    newid = opts.newid
    oo = opts.oo

    o = OO(oo, sizes.mapping)

    if oo:
        seen = o.contigs
        # The leftover contigs not in the oo file
        logging.debug(
            "A total of {0} contigs ({1} in `{2}`)".format(len(sizes), len(seen), oo)
        )

        for ctg, size in sizes.iter_sizes():
            if ctg in seen:
                continue
            o.add(ctg, ctg, size)

    else:
        if newid:
            for ctg, size in sizes.iter_sizes():
                o.add(newid, ctg, size)
        else:
            for scaffold_number, (ctg, size) in enumerate(sizes.iter_sizes()):
                object_id = "scaffold{0:03d}".format(scaffold_number + 1)
                o.add(object_id, ctg, size)

    fw = open(agpfile, "w")
    o.write_AGP(
        fw,
        gapsize=opts.gapsize,
        gaptype=opts.gaptype,
        evidence=opts.evidence,
        phases=phases,
    )
    fw.close()

    joinedfastafile = prefix + ".joined.fasta"
    build([agpfile, fastafile, joinedfastafile])

    return joinedfastafile


def summary(args):
    """
    %prog summary *.fasta

    Report real bases and N's in fastafiles in a tabular report
    """
    from natsort import natsort_key

    p = OptionParser(summary.__doc__)
    p.add_option(
        "--suffix", default="Mb", help="make the base pair counts human readable"
    )
    p.add_option("--ids", help="write the ids that have >= 50% N's")
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    idsfile = opts.ids
    header = "Seqid Real N's Total %_real".split()
    if idsfile:
        idsfile = open(idsfile, "w")
        nids = 0

    data = []
    for fastafile in args:
        for rec in SeqIO.parse(fastafile, "fasta"):
            seqlen = len(rec)
            nns = rec.seq.count("n") + rec.seq.count("N")
            reals = seqlen - nns
            pct = reals * 100.0 / seqlen
            pctreal = "{0:.1f}%".format(pct)
            if idsfile and pct < 50:
                nids += 1
                print(rec.id, file=idsfile)

            data.append((rec.id, reals, nns, seqlen, pctreal))

    data.sort(key=natsort_key)
    ids, reals, nns, seqlen, pctreal = zip(*data)
    reals = sum(reals)
    nns = sum(nns)
    seqlen = sum(seqlen)
    pctreal = "{0:.1f}%".format(reals * 100.0 / seqlen)
    data.append(("Total", reals, nns, seqlen, pctreal))

    write_csv(header, data, sep=" ", filename=opts.outfile, thousands=True)
    if idsfile:
        logging.debug(
            "A total of {0} ids >= 50% N's written to {1}.".format(nids, idsfile.name)
        )
        idsfile.close()

    return reals, nns, seqlen


def format(args):
    """
    %prog format infasta outfasta

    Reformat FASTA file and also clean up names.
    """
    sequential_choices = ("replace", "prefix", "suffix")
    p = OptionParser(format.__doc__)
    p.add_option(
        "--pairs",
        default=False,
        action="store_true",
        help="Add trailing /1 and /2 for interleaved pairs",
    )
    p.add_option(
        "--sequential",
        default=None,
        choices=sequential_choices,
        help="Add sequential IDs",
    )
    p.add_option(
        "--sequentialoffset", default=0, type="int", help="Sequential IDs start at"
    )
    p.add_option(
        "--pad0", default=0, type="int", help="Pad a few zeros in front of sequential"
    )
    p.add_option(
        "--gb",
        default=False,
        action="store_true",
        help="For Genbank ID, get the accession",
    )
    p.add_option("--sep", default=None, help="Split description by certain symbol")
    p.add_option(
        "--index",
        default=0,
        type="int",
        help="Extract i-th field after split with --sep",
    )
    p.add_option(
        "--noversion",
        default=False,
        action="store_true",
        help="Remove the gb trailing version",
    )
    p.add_option("--prefix", help="Prepend prefix to sequence ID")
    p.add_option("--suffix", help="Append suffix to sequence ID")
    p.add_option(
        "--template",
        default=False,
        action="store_true",
        help="Extract `template=aaa dir=x library=m` to `m-aaa/x`",
    )
    p.add_option("--switch", help="Switch ID from two-column file")
    p.add_option(
        "--annotation",
        help="Add functional annotation from two-column file ('ID <--> Annotation')",
    )
    p.add_option("--ids", help="Generate ID conversion table")
    p.add_option(
        "--upper",
        default=False,
        action="store_true",
        help="Convert sequence to upper case",
    )
    p.add_option(
        "--nodesc",
        default=False,
        action="store_true",
        help="Remove description after identifier",
    )
    p.add_option(
        "--minlength", default=0, type="int", help="Minimum sequence length to keep"
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    infasta, outfasta = args
    gb = opts.gb
    pairs = opts.pairs
    prefix = opts.prefix
    suffix = opts.suffix
    noversion = opts.noversion
    sequential = opts.sequential
    sequentialoffset = opts.sequentialoffset
    sep = opts.sep
    idx = opts.index
    mapfile = opts.switch
    annotfile = opts.annotation
    desc = not opts.nodesc
    idsfile = opts.ids
    idsfile = open(idsfile, "w") if idsfile else None
    upper = opts.upper
    minlength = opts.minlength

    if mapfile:
        mapping = DictFile(mapfile, delimiter="\t")
    if annotfile:
        annotation = DictFile(annotfile, delimiter="\t")

    fp = SeqIO.parse(must_open(infasta), "fasta")
    fw = must_open(outfasta, "w")
    nremoved = 0
    for i, rec in enumerate(fp):
        if len(rec) < minlength:
            nremoved += 1
            continue
        origid = rec.id
        description = rec.description.replace(origid, "").strip()
        if sep:
            rec.id = rec.description.split(sep)[idx].strip()
        if gb:
            # gi|262233616|gb|GU123895.1| Coffea arabica clone BAC
            atoms = rec.id.split("|")
            if len(atoms) >= 3:
                rec.id = atoms[3]
            elif len(atoms) == 2:
                rec.id = atoms[1]
        if pairs:
            id = "/1" if (i % 2 == 0) else "/2"
            rec.id += id
        if noversion:
            rec.id = rec.id.rsplit(".", 1)[0]
        if sequential:
            rec.id = "{0:0{1}d}".format(sequentialoffset, opts.pad0)
            if sequential == "prefix":
                rec.id = "{0}-{1}".format(rec.id, origid)
            elif sequential == "suffix":
                rec.id = "{0}-{1}".format(origid, rec.id)
            sequentialoffset += 1
        if opts.template:
            template, dir, lib = [
                x.split("=")[-1] for x in rec.description.split()[1:4]
            ]
            rec.id = "{0}-{1}/{2}".format(lib, template, dir)
        if mapfile:
            if origid in mapping:
                rec.id = mapping[origid]
            else:
                logging.error(
                    "{0} not found in `{1}`. ID unchanged.".format(origid, mapfile)
                )
        if prefix:
            rec.id = prefix + rec.id
        if suffix:
            rec.id += suffix
        if annotfile:
            rec.description = (
                annotation.get(origid, "")
                if not mapfile
                else annotation.get(rec.id, "")
            )
        else:
            rec.description = description if desc else ""
        if idsfile:
            print("\t".join((origid, rec.id)), file=idsfile)
        if upper:
            rec.seq = rec.seq.upper()

        SeqIO.write(rec, fw, "fasta")

    if idsfile:
        logging.debug("Conversion table written to `{0}`.".format(idsfile.name))
        idsfile.close()

    if nremoved:
        logging.debug(
            "Removed {} sequences with length < {}".format(nremoved, minlength)
        )


def print_first_difference(
    arec, brec, ignore_case=False, ignore_N=False, rc=False, report_match=True
):
    """
    Returns the first different nucleotide in two sequence comparisons
    runs both Plus and Minus strand
    """
    plus_match = _print_first_difference(
        arec,
        brec,
        ignore_case=ignore_case,
        ignore_N=ignore_N,
        report_match=report_match,
    )
    if rc and not plus_match:
        logging.debug("trying reverse complement of %s" % brec.id)
        brec.seq = brec.seq.reverse_complement()
        minus_match = _print_first_difference(
            arec,
            brec,
            ignore_case=ignore_case,
            ignore_N=ignore_N,
            report_match=report_match,
        )
        return minus_match

    else:
        return plus_match


def _print_first_difference(
    arec, brec, ignore_case=False, ignore_N=False, report_match=True
):
    """
    Returns the first different nucleotide in two sequence comparisons
    """
    aseq, bseq = arec.seq, brec.seq
    asize, bsize = len(aseq), len(bseq)

    matched = True
    for i, (a, b) in enumerate(zip_longest(aseq, bseq)):
        if ignore_case and None not in (a, b):
            a, b = a.upper(), b.upper()

        if ignore_N and ("N" in (a, b) or "X" in (a, b)):
            continue

        if a != b:
            matched = False
            break

    if i + 1 == asize and matched:
        if report_match:
            printf("[green]Two sequences match")
        match = True
    else:
        printf("[red]Two sequences do not match")

        snippet_size = 20  # show the context of the difference

        printf("[red]Sequence start to differ at position {}:".format(i + 1))

        begin = max(i - snippet_size, 0)
        aend = min(i + snippet_size, asize)
        bend = min(i + snippet_size, bsize)

        printf("[red]{}|{}".format(aseq[begin:i], aseq[i:aend]))
        printf("[red]{}|{}".format(bseq[begin:i], bseq[i:bend]))
        match = False

    return match


def diff(args):
    """
    %prog diff afasta bfasta

    print out whether the records in two fasta files are the same
    """
    from jcvi.utils.table import banner

    p = OptionParser(diff.__doc__)
    p.add_option(
        "--ignore_case",
        default=False,
        action="store_true",
        help="ignore case when comparing sequences",
    )
    p.add_option(
        "--ignore_N",
        default=False,
        action="store_true",
        help="ignore N and X's when comparing sequences",
    )
    p.add_option(
        "--ignore_stop",
        default=False,
        action="store_true",
        help="ignore stop codon when comparing sequences",
    )
    p.add_option(
        "--rc",
        default=False,
        action="store_true",
        help="also consider reverse complement",
    )
    p.add_option(
        "--quiet",
        default=False,
        action="store_true",
        help="don't output comparison details",
    )

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    afasta, bfasta = args

    afastan = len(Fasta(afasta))
    bfastan = len(Fasta(bfasta))

    if afastan == bfastan:
        printf(
            "[green]Two sets contain the same number of sequences ({}, {})".format(
                afastan, bfastan
            )
        )
    else:
        printf(
            "[red]Two sets contain different number of sequences ({}, {})".format(
                afastan, bfastan
            )
        )

    ah = SeqIO.parse(afasta, "fasta")
    bh = SeqIO.parse(bfasta, "fasta")

    problem_ids = []
    for arec, brec in zip(ah, bh):

        if opts.ignore_stop:
            arec.seq = arec.seq.rstrip("*")
            brec.seq = brec.seq.rstrip("*")

        asize, bsize = len(arec), len(brec)

        if not opts.quiet:
            print(banner(str(arec), [str(brec)]))
            if asize == bsize:
                printf("[green]Two sequence size match ({})".format(asize))
            else:
                printf(
                    "[red]Two sequence size do not match ({}, {}})".format(asize, bsize)
                )

        # print out the first place the two sequences diff
        fd = print_first_difference(
            arec,
            brec,
            ignore_case=opts.ignore_case,
            ignore_N=opts.ignore_N,
            rc=opts.rc,
            report_match=not opts.quiet,
        )
        if not fd:
            logging.error("Two sets of sequences differ at `{0}`".format(arec.id))
            problem_ids.append(
                "\t".join(str(x) for x in (arec.id, asize, bsize, abs(asize - bsize)))
            )

    if problem_ids:
        print(red("A total of {0} records mismatch.".format(len(problem_ids))))
        fw = must_open("Problems.ids", "w")
        print("\n".join(problem_ids), file=fw)


def hash_fasta(
    seq, ignore_case=False, ignore_N=False, ignore_stop=False, checksum="MD5"
):
    """
    Generates checksum of input sequence element
    """
    if ignore_stop:
        seq = seq.rstrip("*")
    if ignore_case:
        seq = seq.upper()
    if ignore_N:
        if not all(c.upper() in "ATGCN" for c in seq):
            seq = re.sub("X", "", seq)
        else:
            seq = re.sub("N", "", seq)

    if checksum == "MD5":
        hashed = md5(seq).hexdigest()
    elif checksum == "GCG":
        hashed = seguid(seq)

    return seguid(seq) if checksum == "GCG" else hashlib.sha256(seq)


def identical(args):
    """
    %prog identical *.fasta

    Given multiple fasta files, find all the exactly identical records
    based on the computed md5 hexdigest or GCG checksum of each sequence.

    Output is an N + 1 column file (where N = number of input fasta files).
    If there are duplicates within a given fasta file, they will all be
    listed out in the same row separated by a comma.

    Example output:
    ---------------------------
               tta1.fsa    tta2.fsa
        t0         2131          na
        t1         3420          na
        t2    3836,3847         852
        t3          148         890
        t4          584         614
        t5          623         684
        t6         1281         470
        t7         3367          na
    """
    from jcvi.utils.cbook import AutoVivification

    allowed_checksum = ["MD5", "GCG"]

    p = OptionParser(identical.__doc__)
    p.add_option(
        "--ignore_case",
        default=False,
        action="store_true",
        help="ignore case when comparing sequences",
    )
    p.add_option(
        "--ignore_N",
        default=False,
        action="store_true",
        help="ignore N and X's when comparing sequences",
    )
    p.add_option(
        "--ignore_stop",
        default=False,
        action="store_true",
        help="ignore stop codon when comparing sequences",
    )
    p.add_option(
        "--output_uniq",
        default=False,
        action="store_true",
        help="output uniq sequences in FASTA format",
    )
    p.add_option(
        "--checksum",
        default="MD5",
        choices=allowed_checksum,
        help="specify checksum method",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    d = AutoVivification()
    files = []
    for fastafile in args:
        f = Fasta(fastafile)
        pf = fastafile.rsplit(".", 1)[0]
        files.append(pf)

        logging.debug("Hashing individual elements of {0}".format(fastafile))
        for name, rec in f.iteritems_ordered():
            seq = re.sub(" ", "", str(rec.seq))
            hashed = hash_fasta(
                seq,
                ignore_case=opts.ignore_case,
                ignore_N=opts.ignore_N,
                ignore_stop=opts.ignore_stop,
                checksum=opts.checksum,
            )
            if not d[hashed]:
                d[hashed]["seq"] = seq
                d[hashed]["count"] = 0
            if not d[hashed]["names"][pf]:
                d[hashed]["names"][pf] = set()
            d[hashed]["names"][pf].add(name)

    fw = must_open(opts.outfile, "w")
    if opts.output_uniq:
        uniqfile = "_".join(files) + ".uniq.fasta"
        uniqfw = must_open(uniqfile, "w")

    header = "\t".join(str(x) for x in (args))
    print("\t".join(str(x) for x in ("", header)), file=fw)
    for idx, hashed in enumerate(d.keys()):
        line = []
        line.append("t{0}".format(idx))
        for fastafile in files:
            if fastafile in d[hashed]["names"].keys():
                line.append(",".join(d[hashed]["names"][fastafile]))
                if opts.output_uniq:
                    d[hashed]["count"] += len(d[hashed]["names"][fastafile])
            else:
                line.append("na")
        print("\t".join(line), file=fw)

        if opts.output_uniq:
            seqid = "\t".join(str(x) for x in ("t{0}".format(idx), d[hashed]["count"]))
            rec = SeqRecord(Seq(d[hashed]["seq"]), id=seqid, description="")
            SeqIO.write([rec], uniqfw, "fasta")

    fw.close()
    if opts.output_uniq:
        logging.debug("Uniq sequences written to `{0}`".format(uniqfile))
        uniqfw.close()


QUALSUFFIX = ".qual"


def get_qual(fastafile, suffix=QUALSUFFIX, check=True):
    """
    Check if current folder contains a qual file associated with the fastafile
    """
    qualfile1 = fastafile.rsplit(".", 1)[0] + suffix
    qualfile2 = fastafile + suffix

    if check:
        if op.exists(qualfile1):
            logging.debug("qual file `{0}` found".format(qualfile1))
            return qualfile1
        elif op.exists(qualfile2):
            logging.debug("qual file `{0}` found".format(qualfile2))
            return qualfile2
        else:
            return None

    return qualfile1


def some(args):
    """
    %prog some fastafile listfile outfastafile

    generate a subset of fastafile, based on a list
    """
    from jcvi.utils.cbook import gene_name

    p = OptionParser(some.__doc__)
    p.add_option(
        "--exclude",
        default=False,
        action="store_true",
        help="Output sequences not in the list file",
    )
    p.add_option(
        "--no_strip_names",
        default=False,
        action="store_true",
        help="Do not strip alternative splicing (e.g. At5g06540.1 -> At5g06540)",
    )
    p.add_option(
        "--uniprot", default=False, action="store_true", help="Header is from uniprot"
    )

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(p.print_help())

    strip_names = not opts.no_strip_names
    fastafile, listfile, outfastafile = args
    outfastahandle = must_open(outfastafile, "w")
    qualfile = get_qual(fastafile)

    names = set(open(listfile).read().split())
    if qualfile:
        outqualfile = outfastafile + ".qual"
        outqualhandle = open(outqualfile, "w")
        parser = iter_fasta_qual(fastafile, qualfile)
    else:
        parser = SeqIO.parse(fastafile, "fasta")

    recs = []
    seen = set()
    for rec in parser:
        name = rec.id
        if strip_names:
            name = gene_name(name)

        if name in seen:  # Only report one instance
            continue

        if opts.uniprot:
            name = name.split("|")[-1]

        if opts.exclude:
            if name in names:
                continue
        else:
            if name not in names:
                continue

        recs.append(rec)
        seen.add(name)

    for rec in recs:
        SeqIO.write([rec], outfastahandle, "fasta")
        if qualfile:
            SeqIO.write([rec], outqualhandle, "qual")

    logging.debug("A total of %d records written to `%s`" % (len(recs), outfastafile))


def fastq(args):
    """
    %prog fastq fastafile

    Generate fastqfile by combining fastafile and fastafile.qual.
    Also check --qv option to use a default qv score.
    """
    from jcvi.formats.fastq import FastqLite

    p = OptionParser(fastq.__doc__)
    p.add_option("--qv", type="int", help="Use generic qv value")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    fastqfile = fastafile.rsplit(".", 1)[0] + ".fastq"
    fastqhandle = open(fastqfile, "w")
    num_records = 0

    if opts.qv is not None:
        qv = chr(ord("!") + opts.qv)
        logging.debug("QV char '{0}' ({1})".format(qv, opts.qv))
    else:
        qv = None

    if qv:
        f = Fasta(fastafile, lazy=True)
        for name, rec in f.iteritems_ordered():
            r = FastqLite("@" + name, str(rec.seq).upper(), qv * len(rec.seq))
            print(r, file=fastqhandle)
            num_records += 1

    else:
        qualfile = get_qual(fastafile)
        for rec in iter_fasta_qual(fastafile, qualfile):
            SeqIO.write([rec], fastqhandle, "fastq")
            num_records += 1

    fastqhandle.close()
    logging.debug("A total of %d records written to `%s`" % (num_records, fastqfile))


def pair(args):
    """
    %prog pair fastafile

    Generate .pairs.fasta and .fragments.fasta by matching records
    into the pairs and the rest go to fragments.
    """
    p = OptionParser(pair.__doc__)
    p.set_sep(
        sep=None,
        help="Separator in name to reduce to clone id"
        + "e.g. GFNQ33242/1 use /, BOT01-2453H.b1 use .",
    )
    p.add_option(
        "-m",
        dest="matepairs",
        default=False,
        action="store_true",
        help="generate .matepairs file [often used for Celera Assembler]",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (fastafile,) = args
    qualfile = get_qual(fastafile)

    prefix = fastafile.rsplit(".", 1)[0]
    pairsfile = prefix + ".pairs.fasta"
    fragsfile = prefix + ".frags.fasta"
    pairsfw = open(pairsfile, "w")
    fragsfw = open(fragsfile, "w")

    # TODO: need a class to handle coupled fasta and qual iterating and indexing
    if opts.matepairs:
        matepairsfile = prefix + ".matepairs"
        matepairsfw = open(matepairsfile, "w")

    if qualfile:
        pairsqualfile = pairsfile + ".qual"
        pairsqualhandle = open(pairsqualfile, "w")
        fragsqualfile = fragsfile + ".qual"
        fragsqualhandle = open(fragsqualfile, "w")

    f = Fasta(fastafile)
    if qualfile:
        q = SeqIO.index(qualfile, "qual")

    all_keys = list(f.keys())
    all_keys.sort()
    sep = opts.sep

    if sep:
        key_fun = lambda x: x.split(sep, 1)[0]
    else:
        key_fun = lambda x: x[:-1]

    for key, variants in groupby(all_keys, key=key_fun):
        variants = list(variants)
        paired = len(variants) == 2

        if paired and opts.matepairs:
            print("\t".join(("%s/1" % key, "%s/2" % key)), file=matepairsfw)

        fw = pairsfw if paired else fragsfw
        if qualfile:
            qualfw = pairsqualhandle if paired else fragsqualhandle

        for i, var in enumerate(variants):
            rec = f[var]
            if qualfile:
                recqual = q[var]
            newid = "%s/%d" % (key, i + 1)

            rec.id = newid
            rec.description = ""
            SeqIO.write([rec], fw, "fasta")
            if qualfile:
                recqual.id = newid
                recqual.description = ""
                SeqIO.write([recqual], qualfw, "qual")

    logging.debug("sequences written to `%s` and `%s`" % (pairsfile, fragsfile))
    if opts.matepairs:
        logging.debug("mates written to `%s`" % matepairsfile)


def pairinplace(args):
    """
    %prog pairinplace bulk.fasta

    Pair up the records in bulk.fasta by comparing the names for adjacent
    records. If they match, print to bulk.pairs.fasta, else print to
    bulk.frags.fasta.
    """
    p = OptionParser(pairinplace.__doc__)
    p.add_option(
        "-r",
        dest="rclip",
        default=1,
        type="int",
        help="pair ID is derived from rstrip N chars",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    base = op.basename(fastafile).split(".")[0]

    frags = base + ".frags.fasta"
    pairs = base + ".pairs.fasta"
    if fastafile.endswith(".gz"):
        frags += ".gz"
        pairs += ".gz"

    fragsfw = must_open(frags, "w")
    pairsfw = must_open(pairs, "w")

    N = opts.rclip
    strip_name = lambda x: x[:-N] if N else str

    skipflag = False  # controls the iterator skip
    fastaiter = SeqIO.parse(fastafile, "fasta")
    for a, b in pairwise(fastaiter):

        aid, bid = [strip_name(x) for x in (a.id, b.id)]

        if skipflag:
            skipflag = False
            continue

        if aid == bid:
            SeqIO.write([a, b], pairsfw, "fasta")
            skipflag = True
        else:
            SeqIO.write([a], fragsfw, "fasta")

    # don't forget the last one, when b is None
    if not skipflag:
        SeqIO.write([a], fragsfw, "fasta")

    logging.debug("Reads paired into `%s` and `%s`" % (pairs, frags))


def extract(args):
    """
    %prog extract fasta query

    extract query out of fasta file, query needs to be in the form of
    "seqname", or "seqname:start-stop", or "seqname:start-stop:-"
    """
    p = OptionParser(extract.__doc__)
    p.add_option("--newname", help="Use this new name instead")
    p.add_option(
        "--include",
        default=False,
        action="store_true",
        help="search description line for match",
    )
    p.add_option(
        "--exclude",
        default=False,
        action="store_true",
        help="exclude description that matches",
    )
    p.add_option(
        "--idonly", default=False, action="store_true", help="Only search identifier"
    )
    p.add_option(
        "--bed",
        default=None,
        help="path to bed file to guide extraction by matching seqname",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) == 2:
        fastafile, query = args
    elif len(args) == 1 and opts.bed:
        (fastafile,) = args
        bedaccns = Bed(opts.bed).accns
    else:
        sys.exit(p.print_help())

    if opts.bed:
        fw = must_open(opts.outfile, "w")
        f = Fasta(fastafile)
        for accn in bedaccns:
            try:
                rec = f[accn]
            except:
                logging.error("{0} not found in {1}".format(accn, fastafile))
                continue
            SeqIO.write([rec], fw, "fasta")
        return fw.name

    atoms = query.split(":")
    key = atoms[0]

    assert len(atoms) <= 3, "cannot have more than two ':' in your query"

    pos = ""
    if len(atoms) in (2, 3):
        pos = atoms[1]

    strand = "+"
    if len(atoms) == 3:
        strand = atoms[2]

    assert strand in ("+", "-"), "strand must be either '+' or '-'"

    feature = dict(chr=key)

    if "-" in pos:
        start, stop = pos.split("-")
        try:
            start, stop = int(start), int(stop)
        except ValueError as e:
            logging.error(e)
            sys.exit(p.print_help())

        feature["start"] = start
        feature["stop"] = stop
    else:
        start, stop = None, None

    assert start < stop or None in (
        start,
        stop,
    ), "start must be < stop, you have ({0}, {1})".format(start, stop)
    feature["strand"] = strand

    include, exclude = opts.include, opts.exclude
    # conflicting options, cannot be true at the same time
    assert not (
        include and exclude
    ), "--include and --exclude cannot be on at the same time"
    fw = must_open(opts.outfile, "w")

    if include or exclude:
        f = Fasta(fastafile, lazy=True)
        fi = f.iteritems_ordered if opts.idonly else f.iterdescriptions_ordered
        for k, rec in fi():
            if include and key not in k:
                continue
            if exclude and key in k:
                continue

            seq = Fasta.subseq(rec, start, stop, strand)
            newid = rec.id
            if start is not None:
                newid += ":{0}-{1}:{2}".format(start, stop, strand)

            rec = SeqRecord(seq, id=newid, description=k)
            SeqIO.write([rec], fw, "fasta")
    else:
        f = Fasta(fastafile)
        try:
            seq = f.sequence(feature, asstring=False)
        except AssertionError as e:
            logging.error(e)
            return

        newid = opts.newname or query
        rec = SeqRecord(seq, id=newid, description="")
        SeqIO.write([rec], fw, "fasta")

    return fw.name


def _uniq_rec(fastafile, seq=False):
    """
    Returns unique records
    """
    seen = set()
    for rec in SeqIO.parse(fastafile, "fasta"):
        name = str(rec.seq) if seq else rec.id
        if name in seen:
            logging.debug("ignore {0}".format(rec.id))
            continue
        seen.add(name)
        yield rec


def uniq(args):
    """
    %prog uniq fasta uniq.fasta

    remove fasta records that are the same
    """
    p = OptionParser(uniq.__doc__)
    p.add_option(
        "--seq", default=False, action="store_true", help="Uniqify the sequences"
    )
    p.add_option(
        "-t",
        "--trimname",
        dest="trimname",
        action="store_true",
        default=False,
        help="turn on the defline trim to first space",
    )

    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(p.print_help())

    fastafile, uniqfastafile = args
    fw = must_open(uniqfastafile, "w")
    seq = opts.seq

    for rec in _uniq_rec(fastafile, seq=seq):
        if opts.trimname:
            rec.description = ""
        SeqIO.write([rec], fw, "fasta")


def random(args):
    """
    %prog random fasta 100 > random100.fasta

    Take number of records randomly from fasta
    """
    from random import sample

    p = OptionParser(random.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, N = args
    N = int(N)
    assert N > 0

    f = Fasta(fastafile)
    fw = must_open("stdout", "w")

    for key in sample(f.keys(), N):
        rec = f[key]
        SeqIO.write([rec], fw, "fasta")

    fw.close()


XQUAL = -1000  # default quality for X
NQUAL = 5  # default quality value for N
QUAL = 10  # default quality value
OKQUAL = 15


def modify_qual(rec):
    qv = rec.letter_annotations["phred_quality"]
    for i, (s, q) in enumerate(zip(rec.seq, qv)):
        if s == "X" or s == "x":
            qv[i] = XQUAL
        if s == "N" or s == "x":
            qv[i] = NQUAL
    return rec


def make_qual(fastafile, score=OKQUAL):
    logging.warning("assume qual ({0})".format(score))
    qualfile = fastafile.rsplit(".", 1)[0] + ".qual"
    fw = open(qualfile, "w")
    fasta = Fasta(fastafile, lazy=True)
    score = str(score) + " "
    for entry, size in fasta.itersizes_ordered():
        print(">" + entry, file=fw)
        print(score * size, file=fw)
    fw.close()
    return qualfile


def iter_fasta_qual(fastafile, qualfile, defaultqual=OKQUAL, modify=False):
    """
    used by trim, emits one SeqRecord with quality values in it
    """
    from Bio.SeqIO.QualityIO import PairedFastaQualIterator

    if not qualfile:
        qualfile = make_qual(fastafile, score=defaultqual)

    rec_iter = PairedFastaQualIterator(open(fastafile), open(qualfile))
    for rec in rec_iter:
        yield rec if not modify else modify_qual(rec)


def write_fasta_qual(rec, fastahandle, qualhandle):
    if fastahandle:
        SeqIO.write([rec], fastahandle, "fasta")
    if qualhandle:
        SeqIO.write([rec], qualhandle, "qual")


def trim(args):
    """
    %prog trim fasta.screen newfasta

    take the screen output from `cross_match` (against a vector db, for
    example), then trim the sequences to remove X's. Will also perform quality
    trim if fasta.screen.qual is found. The trimming algorithm is based on
    finding the subarray that maximize the sum
    """

    from jcvi.algorithms.maxsum import max_sum

    p = OptionParser(trim.__doc__)
    p.add_option(
        "-c",
        dest="min_length",
        type="int",
        default=64,
        help="minimum sequence length after trimming",
    )
    p.add_option("-s", dest="score", default=QUAL, help="quality trimming cutoff")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    fastafile, newfastafile = args
    qualfile = get_qual(fastafile)
    newqualfile = get_qual(newfastafile, check=False)

    logging.debug(
        "Trim bad sequence from fasta file `%s` to `%s`" % (fastafile, newfastafile)
    )

    fw = must_open(newfastafile, "w")
    fw_qual = open(newqualfile, "w")

    dropped = trimmed = 0

    for rec in iter_fasta_qual(fastafile, qualfile, modify=True):
        qv = [x - opts.score for x in rec.letter_annotations["phred_quality"]]
        msum, trim_start, trim_end = max_sum(qv)
        score = trim_end - trim_start + 1

        if score < opts.min_length:
            dropped += 1
            continue

        if score < len(rec):
            trimmed += 1
            rec = rec[trim_start : trim_end + 1]

        write_fasta_qual(rec, fw, fw_qual)

    print("A total of %d sequences modified." % trimmed, file=sys.stderr)
    print(
        "A total of %d sequences dropped (length < %d)." % (dropped, opts.min_length),
        file=sys.stderr,
    )

    fw.close()
    fw_qual.close()


def sequin(args):
    """
    %prog sequin inputfasta

    Generate a gapped fasta format with known gap sizes embedded. suitable for
    Sequin submission.

    A gapped sequence represents a newer method for describing non-contiguous
    sequences, but only requires a single sequence identifier. A gap is
    represented by a line that starts with >? and is immediately followed by
    either a length (for gaps of known length) or "unk100" for gaps of unknown
    length. For example, ">?200". The next sequence segment continues on the
    next line, with no separate definition line or identifier. The difference
    between a gapped sequence and a segmented sequence is that the gapped
    sequence uses a single identifier and can specify known length gaps.
    Gapped sequences are preferred over segmented sequences. A sample gapped
    sequence file is shown here:

    >m_gagei [organism=Mansonia gagei] Mansonia gagei NADH dehydrogenase ...
    ATGGAGCATACATATCAATATTCATGGATCATACCGTTTGTGCCACTTCCAATTCCTATTTTAATAGGAA
    TTGGACTCCTACTTTTTCCGACGGCAACAAAAAATCTTCGTCGTATGTGGGCTCTTCCCAATATTTTATT
    >?200
    GGTATAATAACAGTATTATTAGGGGCTACTTTAGCTCTTGC
    TCAAAAAGATATTAAGAGGGGTTTAGCCTATTCTACAATGTCCCAACTGGGTTATATGATGTTAGCTCTA
    >?unk100
    TCAATAAAACTATGGGGTAAAGAAGAACAAAAAATAATTAACAGAAATTTTCGTTTATCTCCTTTATTAA
    TATTAACGATGAATAATAATGAGAAGCCATATAGAATTGGTGATAATGTAAAAAAAGGGGCTCTTATTAC
    """
    p = OptionParser(sequin.__doc__)
    p.add_option("--unk", default=100, type="int", help="The size for unknown gaps")
    p.add_option("--newid", default=None, help="Use this identifier instead")
    p.add_option(
        "--chromosome", default=None, help="Add [chromosome= ] to FASTA header"
    )
    p.add_option("--clone", default=None, help="Add [clone= ] to FASTA header")
    p.set_mingap(default=100)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (inputfasta,) = args
    unk = opts.unk

    outputfasta = inputfasta.rsplit(".", 1)[0] + ".split"
    rec = next(SeqIO.parse(must_open(inputfasta), "fasta"))
    seq = ""
    unknowns, knowns = 0, 0
    for gap, gap_group in groupby(rec.seq, lambda x: x.upper() == "N"):
        subseq = "".join(gap_group)
        if gap:
            gap_length = len(subseq)
            if gap_length == unk:
                subseq = "\n>?unk{0}\n".format(unk)
                unknowns += 1
            elif gap_length >= opts.mingap:
                subseq = "\n>?{0}\n".format(gap_length)
                knowns += 1
        seq += subseq

    fw = must_open(outputfasta, "w")
    id = opts.newid or rec.id
    fastaheader = ">{0}".format(id)
    if opts.chromosome:
        fastaheader += " [chromosome={0}]".format(opts.chromosome)
    if opts.clone:
        fastaheader += " [clone={0}]".format(opts.clone)

    print(fastaheader, file=fw)
    print(seq, file=fw)
    fw.close()
    logging.debug(
        "Sequin FASTA written to `{0}` (gaps: {1} unknowns, {2} knowns).".format(
            outputfasta, unknowns, knowns
        )
    )

    return outputfasta, unknowns + knowns


def remove_small_components(rec, minlen):
    newseq = []
    removed = 0
    for gap, seq in groupby(rec.seq, lambda x: x.upper() == "N"):
        seq = "".join(seq)
        seqlen = len(seq)
        if not gap and seqlen < minlen:
            seq = seqlen * "N"  # Mask small components
            logging.debug("Discard component ({0}) in {1}".format(seqlen, rec.name))
            removed += seqlen
        newseq.append(seq)
    rec.seq = Seq("".join(newseq))
    return removed


def trim_terminal_Ns(rec):
    rec.seq = rec.seq.strip("N")


def normalize_gaps(rec, gapsize):
    newseq = []
    normalized = 0
    NN = gapsize * "N"
    for gap, seq in groupby(rec.seq, lambda x: x.upper() == "N"):
        seq = "".join(seq)
        if gap:
            seq = NN
            normalized += 1
        newseq.append(seq)
    rec.seq = Seq("".join(newseq))
    return normalized


def tidy(args):
    """
    %prog tidy fastafile

    Trim terminal Ns, normalize gap sizes and remove small components.
    """
    p = OptionParser(tidy.__doc__)
    p.add_option(
        "--gapsize",
        dest="gapsize",
        default=0,
        type="int",
        help="Set all gaps to the same size",
    )
    p.add_option(
        "--minlen",
        dest="minlen",
        default=100,
        type="int",
        help="Minimum component size",
    )

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    gapsize = opts.gapsize
    minlen = opts.minlen

    tidyfastafile = fastafile.rsplit(".", 1)[0] + ".tidy.fasta"
    fw = must_open(tidyfastafile, "w")

    removed = normalized = 0
    fasta = Fasta(fastafile, lazy=True)
    for name, rec in fasta.iteritems_ordered():
        rec.seq = rec.seq.upper()
        if minlen:
            removed += remove_small_components(rec, minlen)
        trim_terminal_Ns(rec)
        if gapsize:
            normalized += normalize_gaps(rec, gapsize)

        if len(rec) == 0:
            logging.debug("Drop seq {0}".format(rec.id))
            continue
        SeqIO.write([rec], fw, "fasta")

    # Print statistics
    if removed:
        logging.debug("Total discarded bases: {0}".format(removed))
    if normalized:
        logging.debug("Gaps normalized: {0}".format(normalized))

    logging.debug("Tidy FASTA written to `{0}`.".format(tidyfastafile))
    fw.close()

    return tidyfastafile


def write_gaps_worker(rec):
    start = 0
    seq = rec.seq.upper()
    output = []
    for gap, seq in groupby(seq, lambda x: x == "N"):
        seq = "".join(seq)
        current_length = len(seq)
        object_beg = start + 1
        object_end = start + current_length
        if gap:
            s = "\t".join(str(x) for x in (rec.id, object_beg - 1, object_end))
            output.append(s)
        start += current_length

    return "\n".join(output)


def write_gaps_bed(inputfasta, prefix, mingap, cpus):
    from jcvi.apps.grid import WriteJobs
    from jcvi.formats.bed import sort

    bedfile = prefix + ".gaps.bed"
    f = Fasta(inputfasta)
    recs = list(rec for k, rec in f.iteritems())
    pool = WriteJobs(write_gaps_worker, recs, bedfile, cpus=cpus)
    pool.run()

    sort([bedfile, "-i"])

    bed = Bed(bedfile)
    nbedfile = prefix + ".{0}N.bed".format(mingap)

    gapnum = 0
    fw = open(nbedfile, "w")
    for b in bed:
        if b.span < mingap:
            continue
        gapnum += 1
        gapname = "gap.{0:05d}".format(gapnum)
        print("\t".join(str(x) for x in (b, gapname, b.span)), file=fw)

    shutil.move(nbedfile, bedfile)
    logging.debug("Write gap (>={0}bp) locations to `{1}`.".format(mingap, bedfile))


def gaps(args):
    """
    %prog gaps fastafile

    Print out a list of gaps in BED format (.gaps.bed).
    """
    from jcvi.formats.sizes import agp
    from jcvi.formats.agp import mask, build

    p = OptionParser(gaps.__doc__)
    p.add_option(
        "--split", default=False, action="store_true", help="Generate .split.fasta"
    )
    p.set_mingap(default=100)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (inputfasta,) = args
    mingap = opts.mingap
    split = opts.split
    prefix = inputfasta.rsplit(".", 1)[0]
    bedfile = prefix + ".gaps.bed"

    if need_update(inputfasta, bedfile):
        write_gaps_bed(inputfasta, prefix, mingap, opts.cpus)

    if split:
        splitfile = prefix + ".split.fasta"
        oagpfile = prefix + ".splitobject.agp"
        cagpfile = prefix + ".splitcomponent.agp"

        if need_update((inputfasta, bedfile), splitfile):

            sizesagpfile = agp([inputfasta])

            maskedagpfile = mask([sizesagpfile, bedfile, "--splitobject"])
            shutil.move(maskedagpfile, oagpfile)
            logging.debug("AGP file written to `{0}`.".format(oagpfile))

            maskedagpfile = mask([sizesagpfile, bedfile, "--splitcomponent"])
            shutil.move(maskedagpfile, cagpfile)
            logging.debug("AGP file written to `{0}`.".format(cagpfile))

            build([oagpfile, inputfasta, splitfile])
            os.remove(sizesagpfile)

        return splitfile, oagpfile, cagpfile


if __name__ == "__main__":
    main()
