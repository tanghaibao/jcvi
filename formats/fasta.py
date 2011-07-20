"""
Wrapper for biopython Fasta, add option to parse sequence headers
"""

import sys
import os.path as op
import logging

from random import sample
from optparse import OptionParser
from itertools import groupby, izip_longest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcvi.formats.base import BaseFile, DictFile, must_open
from jcvi.utils.cbook import human_size
from jcvi.utils.table import banner
from jcvi.apps.base import ActionDispatcher, debug
from jcvi.apps.console import print_red, print_green
debug()


class Fasta (BaseFile, dict):

    def __init__(self, filename, index=True, key_function=None, lazy=False):
        super(Fasta, self).__init__(filename)
        self.key_function = key_function

        if lazy:  # do not incur the overhead
            return

        if index:
            self.index = SeqIO.index(filename, "fasta",
                    key_function=key_function)
        else:
            # SeqIO.to_dict expects a different key_function that operates on
            # the SeqRecord instead of the raw string
            _key_function = (lambda rec: key_function(rec.description)) if \
                    key_function else None
            self.index = SeqIO.to_dict(SeqIO.parse(must_open(filename), "fasta"),
                    key_function=_key_function)

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
        for k in self.index.iterkeys():
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

    def iterkeys_ordered(self):
        for k, rec in self.iteritems_ordered():
            yield k

    def itersizes_ordered(self):
        for k, rec in self.iteritems_ordered():
            yield k, len(rec)

    def iter_annotations(self):
        for k, rec in self.iteritems_ordered():
            description = rec.description
            yield k, description

    @classmethod
    def subseq(cls, fasta, start=None, stop=None, strand=None):
        """
        Take Bio.SeqRecord and slice "start:stop" from it, does proper
        index and error handling
        """
        start = start - 1 if start is not None else 0
        stop = stop if stop is not None else len(fasta)

        assert start >= 0, "start (%d) must >= 0" % (start + 1)

        assert stop <= len(fasta), \
                ("stop (%d) must be <= " + \
                "length of `%s` (%d)") % (stop, fasta.id, len(fasta))

        seq = fasta.seq[start:stop]

        if strand in (-1, '-1', '-'):
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

        assert 'chr' in f, "`chr` field required"
        name = f['chr']

        assert name in self, "feature: %s not in `%s`" % \
                (f, self.filename)

        fasta = self[f['chr']]

        seq = Fasta.subseq(fasta,
                f.get('start'), f.get('stop'), f.get('strand'))

        if asstring:
            return str(seq)

        return seq


def main():

    actions = (
        ('extract', 'given fasta file and seq id, retrieve the sequence ' + \
                    'in fasta format'),
        ('summary', "report the real no of bases and N's in fastafiles"),
        ('uniq', 'remove records that are the same'),
        ('ids', 'generate a list of header without the >'),
        ('format', 'trim accession id to the first space or switch id ' + \
                   'based on 2-column mapping file'),
        ('pool', 'pool a bunch of fastafiles together and add prefix'),
        ('random', 'randomly take some records'),
        ('diff', 'check if two fasta records contain same information'),
        ('trim', 'given a cross_match screened fasta, trim the sequence'),
        ('sort', 'sort the records by IDs, sizes, etc.'),
        ('pair', 'sort paired reads to .pairs, rest to .fragments'),
        ('fastq', 'combine fasta and qual to create fastq file'),
        ('tidy', 'normalize gap sizes and remove small components in fasta'),
        ('sequin', 'generate a gapped fasta file for sequin submission'),
        ('gaps', 'print out a list of gap sizes within sequences'),
        ('join', 'concatenate a list of seqs and add gaps in between'),
        ('some', 'include or exclude a list of records (also performs on ' + \
                 '.qual file if available)'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def pool(args):
    """
    %prog pool fastafiles

    Pool a bunch of FASTA files, and add prefix to each record based on
    filenames.
    """
    p = OptionParser(pool.__doc__)

    if len(args) < 1:
        sys.exit(not p.print_help())

    for fastafile in args:
        pf = op.basename(fastafile).split(".")[0].split("_")[0]
        prefixopt = "--prefix={0}_".format(pf)
        format([fastafile, "stdout", prefixopt])


def ids(args):
    """
    %prog ids fastafiles

    Generate the FASTA headers without the '>'.
    """
    p = OptionParser(ids.__doc__)

    if len(args) < 1:
        sys.exit(not p.print_help())

    for fastafile in args:
        f = Fasta(fastafile, lazy=True)
        for key in f.iterkeys_ordered():
            print key


def sort(args):
    """
    %prog sort fastafile

    Sort a list of sequences and output with sorted IDs, etc.
    """
    p = OptionParser(sort.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    fastafile, = args
    sortedfastafile = fastafile.rsplit(".", 1)[0] + ".sorted.fasta"

    f = Fasta(fastafile)
    fw = must_open(sortedfastafile, "w")
    for key in sorted(f.iterkeys()):
        rec = f[key]
        SeqIO.write([rec], fw, "fasta")

    logging.debug("Sorted file written to `{0}`.".format(sortedfastafile))
    fw.close()


def join(args):
    """
    %prog join fastafile

    Concatenate a list of seqs and add gaps in between
    """
    p = OptionParser(join.__doc__)
    p.add_option("--minctgsize", default=0, type="int",
            help="minimum contig size allowed [default: %default]")
    p.add_option("--gapsize", dest="gapsize", default=100, type="int",
            help="number of N's in between the sequences [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args

    gap = opts.gapsize * 'N'
    seq = gap.join(str(x.seq) for x in SeqIO.parse(fastafile, "fasta") \
            if len(x.seq) >= opts.minctgsize)
    rec = SeqRecord(Seq(seq), id="chr0", description="")

    fw = sys.stdout
    SeqIO.write([rec], fw, "fasta")


def summary(args):
    """
    %prog summary *.fasta

    Report real bases and N's in fastafiles in a tabular report
    """
    p = OptionParser(summary.__doc__)
    p.add_option("--suffix", dest="suffix", default="Mb",
            help="make the base pair counts human readable [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(p.print_help())

    header = "Seqid|Real|N's|Total|% real"
    print header.replace('|', '\t')
    _human_size = lambda x: human_size(x, precision=2,
            target=opts.suffix)

    total_nns = 0
    total_reals = 0

    data = []
    for fastafile in args:
        for rec in SeqIO.parse(fastafile, "fasta"):
            seqlen = len(rec)
            nns = rec.seq.count('n') + rec.seq.count('N')
            reals = seqlen - nns

            data.append((rec.id, reals, nns, seqlen))

    idsx, total_reals, total_nns, total_seqlen = zip(*data)
    total_reals = sum(total_reals)
    total_nns = sum(total_nns)
    total_seqlen = sum(total_seqlen)
    data.append(("Total", total_reals, total_nns, total_seqlen))

    for id, reals, nns, seqlen in data:
        pctreal = reals * 100. / seqlen
        seqlen = _human_size(seqlen)
        nns = _human_size(nns)
        reals = _human_size(reals)

        print "{0}\t{1}\t{2}\t{3}\t{4:.2f}%".format(id, reals,
                nns, seqlen, pctreal)


def format(args):
    """
    %prog format infasta outfasta

    Reformat FASTA file and also clean up names
    """
    p = OptionParser(format.__doc__)
    p.add_option("--pairs", dest="pairs", default=False, action="store_true",
            help="If input reads are pairs, add trailing /1 and /2 "
            "[default: %default]")
    p.add_option("--sequential", default=False, action="store_true",
            help="Add sequential IDs [default: %default]")
    p.add_option("--pad0", default=6, type="int",
            help="Pad a few zeros in front of sequential [default: %default]")
    p.add_option("--gb", dest="gb", default=False, action="store_true",
            help="For Genbank ID, get the accession [default: %default]")
    p.add_option("--until", default=None,
            help="Get the names until certain symbol [default: %default]")
    p.add_option("--noversion", dest="noversion", default=False,
            action="store_true", help="remove the gb trailing version "
            "[default: %default]")
    p.add_option("--prefix", dest="prefix", default="",
            help="Prepend prefix to the sequence ID [default: '%default']")
    p.add_option("--index", dest="index", default=0, type="int",
            help="Extract i-th field in the description [default: %default]")
    p.add_option("--switch", dest="switch", default=None,
            help="Switch sequence ID based on 2-column mapping file [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    infasta, outfasta = args
    gb = opts.gb
    until = opts.until
    pairs = opts.pairs
    prefix = opts.prefix
    noversion = opts.noversion
    sequential = opts.sequential
    idx = opts.index
    mapfile = opts.switch

    if mapfile:
        mapping = DictFile(mapfile, delimiter="\t")

    fw = must_open(outfasta, "w")
    for i, rec in enumerate(SeqIO.parse(infasta, "fasta")):
        if until:
            rec.id = rec.id.split(until, 1)[0]
        if gb:
            # gi|262233616|gb|GU123895.1| Coffea arabica clone BAC
            atoms = rec.id.split("|")
            if len(atoms) >= 4:
                rec.id = atoms[3]
            elif len(atoms) == 2:
                rec.id = atoms[1]
        if pairs:
            id = "/1" if (i % 2 == 0) else "/2"
            rec.id += id
        if noversion:
            rec.id = rec.id.rsplit(".", 1)[0]
        if sequential:
            rec.id = "{0:0{1}d}".format(i + 1, opts.pad0)
        if prefix:
            rec.id = prefix + rec.id
        if idx:
            rec.id = rec.description.split()[idx]
        if mapfile:
            if rec.id in mapping:
                rec.id = mapping[rec.id]
            else:
                logging.error("{0} not found in `{1}`. ID unchanged.".\
                        format(rec.id, mapfile))
        rec.description = ""

        SeqIO.write(rec, fw, "fasta")


def get_first_rec(fastafile):
    """
    Returns the first record in the fastafile
    """
    f = list(SeqIO.parse(fastafile, "fasta"))

    if len(f) > 1:
        logging.debug("{0} records found in {1}, using the first one".\
                format(len(f), fastafile))

    return f[0]


def print_first_difference(arec, brec, ignore_case=False, ignore_N=False,
        rc=False):
    """
    Returns the first different nucleotide in two sequence comparisons
    runs both Plus and Minus strand
    """
    plus_match = _print_first_difference(arec, brec, ignore_case=ignore_case,
            ignore_N=ignore_N)
    if rc:
        logging.debug("trying reverse complement of %s" % brec.id)
        brec.seq = brec.seq.reverse_complement()
        minus_match = _print_first_difference(arec, brec,
                ignore_case=ignore_case, ignore_N=ignore_N)
        return any((plus_match, minus_match))

    else:
        return plus_match


def _print_first_difference(arec, brec, ignore_case=False, ignore_N=False):
    """
    Returns the first different nucleotide in two sequence comparisons
    """
    aseq, bseq = arec.seq, brec.seq
    asize, bsize = len(aseq), len(bseq)

    for i, (a, b) in enumerate(izip_longest(aseq, bseq)):
        if ignore_case and None not in (a, b):
            a, b = a.upper(), b.upper()

        if (ignore_N and 'N' in (a, b)):
            continue
        if a != b:
            break

    if i + 1 == asize and i + 1 == bsize:
        print_green("Two sequences match")
        match = True
    else:
        print_red("Two sequences do not match")

        snippet_size = 20  # show the context of the difference

        print_red("Sequence start to differ at position %d:" % (i + 1))

        begin = max(i - snippet_size, 0)
        aend = min(i + snippet_size, asize)
        bend = min(i + snippet_size, bsize)

        print_red(aseq[begin:i] + "|" + aseq[i:aend])
        print_red(bseq[begin:i] + "|" + bseq[i:bend])
        match = False

    return match


def diff(args):
    """
    %prog diff afasta bfasta

    print out whether the records in two fasta files are the same
    """
    p = OptionParser(diff.__doc__)
    p.add_option("--ignore_case", dest="ignore_case",
            default=False, action="store_true",
            help="ignore case when comparing sequences")
    p.add_option("--ignore_N", dest="ignore_N",
            default=False, action="store_true",
            help="ignore N's when comparing sequences")
    p.add_option("--rc", dest="rc",
            default=False, action="store_true",
            help="also consider reverse complement")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    afasta, bfasta = args
    #arec = get_first_rec(afasta)
    #brec = get_first_rec(bfasta)

    afastan = len(Fasta(afasta))
    bfastan = len(Fasta(bfasta))

    if afastan == bfastan:
        print_green("Two sets contain the same number of sequences ({0}, {1})".\
                format(afastan, bfastan))
    else:
        print_red("Two sets contain different number of sequences ({0}, {1})".\
                format(afastan, bfastan))

    ah = SeqIO.parse(afasta, "fasta")
    bh = SeqIO.parse(bfasta, "fasta")

    for arec, brec in zip(ah, bh):
        print banner((arec, brec))
        asize, bsize = len(arec), len(brec)

        if asize == bsize:
            print_green("Two sequence size match (%d)" % asize)
        else:
            print_red("Two sequence size do not match (%d, %d)" % (asize, bsize))

        # print out the first place the two sequences diff
        fd = print_first_difference(arec, brec, ignore_case=opts.ignore_case,
                ignore_N=opts.ignore_N, rc=opts.rc)
        if not fd:
            logging.error("Two sets of sequences differ at `{0}`".format(arec.id))
            break


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
            qualfile = qualfile1
        elif op.exists(qualfile2):
            logging.debug("qual file `{0}` found".format(qualfile2))
            qualfile = qualfile2
        else:
            logging.warning("qual file not found")
            return None

    return qualfile1


def some(args):
    """
    %prog some fastafile listfile outfastafile

    generate a subset of fastafile, based on a list
    """
    p = OptionParser(some.__doc__)
    p.add_option("--exclude", dest="exclude",
            default=False, action="store_true",
            help="output sequences not in the list file")

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(p.print_help())

    fastafile, listfile, outfastafile = args
    outfastahandle = must_open(outfastafile, "w")
    qualfile = get_qual(fastafile)

    names = set(x.strip() for x in open(listfile))
    if qualfile:
        outqualfile = outfastafile + ".qual"
        outqualhandle = open(outqualfile, "w")

    parser = iter_fasta_qual(fastafile, qualfile)

    num_records = 0
    for rec in parser:
        name = rec.id
        if opts.exclude:
            if name in names:
                continue
        else:
            if name not in names:
                continue

        rec.description = ""
        SeqIO.write([rec], outfastahandle, "fasta")
        if qualfile:
            SeqIO.write([rec], outqualhandle, "qual")

        num_records += 1

    logging.debug("A total of %d records written to `%s`" % \
            (num_records, outfastafile))


def fastq(args):
    """
    %prog fastq fastafile

    generate fastqfile by combining fastafile and fastafile.qual
    """
    p = OptionParser(fastq.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    fastafile, = args
    qualfile = get_qual(fastafile)
    fastqfile = fastafile.rsplit(".", 1)[0] + ".fastq"
    fastqhandle = open(fastqfile, "w")

    num_records = 0
    for rec in iter_fasta_qual(fastafile, qualfile):
        SeqIO.write([rec], fastqhandle, "fastq")
        num_records += 1
    fastqhandle.close()

    logging.debug("A total of %d records written to `%s`" % \
            (num_records, fastqfile))


def pair(args):
    """
    %prog pair fastafile

    Generate .pairs.fasta and .fragments.fasta by matching records with /1/2
    into the pairs and the rest go to fragments.
    """
    p = OptionParser(pair.__doc__)
    p.add_option("-d", dest="separator", default=None,
            help="separater in the name field to reduce to the same clone " +\
                 "[e.g. GFNQ33242/1 use /, BOT01-2453H.b1 use .]" +\
                 "[default: trim until last char]")
    p.add_option("-m", dest="matepairs", default=False, action="store_true",
            help="generate .matepairs file [often used for Celera Assembler]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    fastafile = args[0]
    qualfile = get_qual(fastafile)

    prefix = fastafile.rsplit(".", 1)[0]
    pairsfile = prefix + ".pairs.fasta"
    fragsfile = prefix + ".frags.fasta"
    pairsfw = open(pairsfile, "w")
    fragsfw = open(fragsfile, "w")

    #TODO: need a class to handle coupled fasta and qual iterating and indexing
    if opts.matepairs:
        matepairsfile = prefix + ".matepairs"
        matepairsfw = open(matepairsfile, "w")

    if qualfile:
        pairsqualfile = pairsfile + ".qual"
        pairsqualhandle = open(pairsqualfile, "w")
        fragsqualfile = fragsfile + ".qual"
        fragsqualhandle = open(fragsqualfile, "w")

    f = Fasta(args[0])
    if qualfile:
        q = SeqIO.index(qualfile, "qual")

    all_keys = list(f.iterkeys())
    all_keys.sort()
    sep = opts.separator

    if sep:
        key_fun = lambda x: x.split(sep, 1)[0]
    else:
        key_fun = lambda x: x[:-1]

    for key, variants in groupby(all_keys, key=key_fun):
        variants = list(variants)
        paired = (len(variants) == 2)

        if paired and opts.matepairs:
            print >> matepairsfw, "\t".join(("%s/1" % key, "%s/2" % key))

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

    logging.debug("sequences written to `%s` and `%s`" % \
            (pairsfile, fragsfile))
    if opts.matepairs:
        logging.debug("mates written to `%s`" % matepairsfile)


def extract(args):
    """
    %prog extract fasta query

    extract query out of fasta file, query needs to be in the form of
    "seqname", or "seqname:start-stop", or "seqname:start-stop:-"
    """
    p = OptionParser(extract.__doc__)
    p.add_option('--include', dest="include",
            default=False, action="store_true",
            help="search description line for match, use 'all' for all " +\
            "records [default: %default]")
    p.add_option('--exclude', dest="exclude",
            default=False, action="store_true",
            help="exclude description that matches [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    fastafile, query = args

    atoms = query.split(":")
    key = atoms[0]

    assert len(atoms) <= 3, "cannot have more than two ':' in your query"

    pos = ""
    if len(atoms) in (2, 3):
        pos = atoms[1]

    strand = "+"
    if len(atoms) == 3:
        strand = atoms[2]

    assert strand in ('+', '-'), "strand must be either '+' or '-'"

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

    assert start < stop or None in (start, stop), \
            "start must be < stop, you have ({0}, {1})".format(start, stop)
    feature["strand"] = strand

    include, exclude = opts.include, opts.exclude
    # conflicting options, cannot be true at the same time
    assert not (include and exclude), "--include and --exclude cannot be "\
            "on at the same time"

    f = Fasta(fastafile)

    if include or exclude:
        for k in f.keys():
            if include and key != "all" and key not in k:
                continue
            if exclude and key in k:
                continue

            rec = f[k]
            seq = Fasta.subseq(rec, start, stop, strand)
            newid = rec.id
            if start is not None:
                newid += ":{0}-{1}:{2}".format(start, stop, strand)

            rec = SeqRecord(seq, id=newid, description="")
            SeqIO.write([rec], sys.stdout, "fasta")
    else:
        try:
            seq = f.sequence(feature, asstring=False)
        except AssertionError as e:
            logging.error(e)
            return

        rec = SeqRecord(seq, id=query, description="")
        SeqIO.write([rec], sys.stdout, "fasta")


def _uniq_rec(fastafile):
    """
    Returns unique records
    """
    seen = set()
    for rec in SeqIO.parse(fastafile, "fasta"):
        name = rec.id
        if name in seen:
            logging.debug("ignore %s" % name)
            continue
        seen.add(name)
        yield rec


def uniq(args):
    """
    %prog uniq fasta uniq.fasta

    remove fasta records that are the same
    """
    p = OptionParser(uniq.__doc__)
    p.add_option("-t", "--trimname", dest="trimname",
            action="store_true", default=False,
            help="turn on the defline trim to first space [default: %default]")

    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(p.print_help())

    fastafile, uniqfastafile = args
    fw = must_open(uniqfastafile, "w")

    for rec in _uniq_rec(fastafile):
        if opts.trimname:
            rec.description = ""
        SeqIO.write([rec], fw, "fasta")


def random(args):
    """
    %prog random fasta 100 > random100.fasta

    Take number of records randomly from fasta
    """
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
    qv = rec.letter_annotations['phred_quality']
    for i, (s, q) in enumerate(zip(rec.seq, qv)):
        if s == 'X' or s == 'x':
            qv[i] = XQUAL
        if s == 'N' or s == 'x':
            qv[i] = NQUAL
    return rec


def iter_fasta_qual(fastafile, qualfile, defaultqual=OKQUAL, modify=False):
    """
    used by trim, emits one SeqRecord with quality values in it
    """
    fastahandle = SeqIO.parse(fastafile, "fasta")

    if qualfile:
        qualityhandle = SeqIO.parse(qualfile, "qual")
        for rec, rec_qual in zip(fastahandle, qualityhandle):
            assert len(rec) == len(rec_qual)
            rec.letter_annotations['phred_quality'] = \
                rec_qual.letter_annotations['phred_quality']
            yield rec if not modify else modify_qual(rec)

    else:
        logging.warning("assume qual ({0})".format(defaultqual))
        for rec in fastahandle:
            rec.letter_annotations['phred_quality'] = [defaultqual] * len(rec)
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
    p.add_option("-c", dest="min_length", type="int", default=64,
            help="minimum sequence length after trimming")
    p.add_option("-s", dest="score", default=QUAL,
            help="quality trimming cutoff [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    fastafile, newfastafile = args
    qualfile = get_qual(fastafile)
    newqualfile = get_qual(newfastafile, check=False)

    logging.debug("Trim bad sequence from fasta file `%s` to `%s`" % \
            (fastafile, newfastafile))

    fw = must_open(newfastafile, "w")
    fw_qual = open(newqualfile, "w")

    dropped = trimmed = 0

    for rec in iter_fasta_qual(fastafile, qualfile, modify=True):
        qv = [x - opts.score for x in \
                rec.letter_annotations["phred_quality"]]
        msum, trim_start, trim_end = max_sum(qv)
        score = trim_end - trim_start + 1

        if score < opts.min_length:
            dropped += 1
            continue

        if score < len(rec):
            trimmed += 1
            rec = rec[trim_start:trim_end + 1]

        write_fasta_qual(rec, fw, fw_qual)

    print >>sys.stderr, "A total of %d sequences modified." % trimmed
    print >>sys.stderr, "A total of %d sequences dropped (length < %d)." % \
        (dropped, opts.min_length)

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
    p.add_option("--mingap", dest="mingap", default=100, type="int",
            help="The minimum size of a gap to split [default: %default]")
    p.add_option("--unk", default=100, type="int",
            help="The size for unknown gaps [default: %default]")
    p.add_option("--newid", default=None,
            help="Use this identifier instead [default: %default]")
    p.add_option("--chromosome", default=None,
            help="Add [chromosome= ] to FASTA header [default: %default]")
    p.add_option("--clone", default=None,
            help="Add [clone= ] to FASTA header [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    inputfasta, = args
    unk = opts.unk

    outputfasta = inputfasta.rsplit(".", 1)[0] + ".split"
    rec = SeqIO.parse(must_open(inputfasta), "fasta").next()
    seq = ""
    unknowns, knowns = 0, 0
    for gap, gap_group in groupby(rec.seq, lambda x: x.upper() == 'N'):
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

    print >> fw, fastaheader
    print >> fw, seq
    fw.close()
    logging.debug("Sequin FASTA written to `{0}` (gaps: {1} unknowns, {2} knowns).".\
            format(outputfasta, unknowns, knowns))

    return outputfasta, unknowns + knowns


def tidy(args):
    """
    %prog tidy fastafile

    Normalize gap sizes (default 100 N's) and remove small components (less than
    100 nucleotides).
    """
    p = OptionParser(tidy.__doc__)
    p.add_option("--gapsize", dest="gapsize", default=100, type="int",
            help="Set all gaps to the same size [default: %default]")
    p.add_option("--minlen", dest="minlen", default=100, type="int",
            help="Minimum component size [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    gapsize = opts.gapsize
    minlen = opts.minlen
    tidyfastafile = fastafile.rsplit(".", 1)[0] + ".tidy.fasta"
    fw = must_open(tidyfastafile, "w")
    normalized_gap = "N" * gapsize

    for rec in SeqIO.parse(fastafile, "fasta"):
        newseq = ""
        dangle_gaps = 0
        for gap, seq in groupby(rec.seq, lambda x: x.upper() == 'N'):
            seq = "".join(seq)
            seqlen = len(seq)
            msg = None
            if gap:
                nsize = max(gapsize - dangle_gaps, 0)
                if seqlen < 10:
                    if nsize > seqlen:
                        nsize = seqlen
                    dangle_gaps += seqlen
                else:
                    if seqlen != gapsize:
                        msg = "Normalize gap size ({0}) to {1}" \
                                .format(seqlen, nsize)
                    dangle_gaps = gapsize

                newseq += nsize * 'N'
            else:
                if seqlen < minlen:
                    msg = "Discard component ({0})".format(seqlen)
                else:
                    newseq += seq
                    # Discarding components might cause flank gaps to merge
                    # should be handled in dangle_gaps, which is only reset when
                    # seeing an actual sequence
                    dangle_gaps = 0

            if msg:
                msg = rec.id + ": " + msg
                logging.info(msg)

        newseq = newseq.strip('N')
        ngaps = newseq.count(normalized_gap)
        if ngaps == 0:
            logging.debug("{0}: is now a Phase 3 sequence.".format(rec.id))
            print "\t".join((rec.id, "3"))

        rec.seq = Seq(newseq)

        SeqIO.write([rec], fw, "fasta")


def gaps(args):
    """
    %prog gaps fastafile

    Print out a list of gaps per sequence record
    """
    p = OptionParser(gaps.__doc__)
    p.add_option("--mingap", dest="mingap", default=10, type="int",
            help="The minimum size of a gap to split [default: %default]")
    p.add_option("--agp", dest="agp", default=False, action="store_true",
            help="Generate AGP file to show components [default: %default]")
    p.add_option("--split", dest="split", default=False, action="store_true",
            help="Generate .split.fasta for the non-gap sequences "
            "[default: %default]")
    p.add_option("--bed", dest="bed", default=False, action="store_true",
            help="Generate .gaps.bed with gap positions [default: %default]")
    p.add_option("--log", dest="log", default=False, action="store_true",
            help="Generate .log with detailed gap positions [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    inputfasta, = args
    mingap = opts.mingap
    agp = opts.agp
    bed = opts.bed
    log = opts.log
    prefix = inputfasta.rsplit(".", 1)[0]
    logfile = prefix + ".log"
    if agp:
        agpfile = prefix + ".agp"
        fwagp = open(agpfile, "w")
        logging.debug("Write AGP to `{0}`.".format(agpfile))
    if opts.split:
        splitfile = prefix + ".split.fasta"
        fwsplit = open(splitfile, "w")
        logging.debug("Write splitted FASTA to `{0}`.".format(splitfile))
    if bed:
        bedfile = prefix + ".gaps.bed"
        fwbed = open(bedfile, "w")
        logging.debug("Write gap locations to `{0}`.".format(bedfile))
    if log:
        fwlog = open(logfile, "w")
        logging.debug("Write gap locations to `{0}`.".format(logfile))

    for rec in SeqIO.parse(inputfasta, "fasta"):
        allgaps = []
        start = 0
        object = rec.id
        component_number = part_number = 0
        for gap, seq in groupby(rec.seq, lambda x: x.upper() == 'N'):
            seq = "".join(seq)
            current_length = len(seq)
            object_beg = start + 1
            object_end = start + current_length
            part_number += 1
            if gap:
                if current_length >= opts.mingap:
                    allgaps.append((current_length, start))
                if agp:
                    component_type = "N"
                    gap_length = current_length
                    gap_type = "fragment"
                    linkage = "yes"
                    empty = ""
                    print >> fwagp, "\t".join(str(x) for x in (object, object_beg,
                        object_end, part_number, component_type, gap_length,
                        gap_type, linkage, empty))
                if bed and len(seq) >= mingap:
                    print >> fwbed, "\t".join(str(x) for x in (object,
                        object_beg, object_end, "gap"))

            else:
                component_id = "{0}_{1}".format(object, component_number)
                component_number += 1
                if agp:
                    component_type = "W"
                    component_beg = 1
                    component_end = current_length
                    orientation = "+"
                    print >> fwagp, "\t".join(str(x) for x in (object, object_beg,
                        object_end, part_number, component_type, component_id,
                        component_beg, component_end, orientation))
                if opts.split:
                    splitrec = SeqRecord(Seq(seq), id=component_id,
                            description="")
                    SeqIO.write([splitrec], fwsplit, "fasta")

            start += current_length

        if allgaps:
            lengths, starts = zip(*allgaps)
            gap_description = ",".join(str(x) for x in lengths)
            starts = ",".join(str(x) for x in starts)
        else:
            gap_description = starts = "no gaps"

        if log:
            print >> fwlog, "\t".join((rec.id, gap_description, starts))


if __name__ == '__main__':
    main()
