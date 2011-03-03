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
from Bio.SeqRecord import SeqRecord

from jcvi.formats.base import BaseFile
from jcvi.utils.cbook import human_size
from jcvi.apps.base import ActionDispatcher, debug
from jcvi.apps.console import tabular, print_red, print_green
debug()


class Fasta (BaseFile, dict):

    def __init__(self, filename, index=True, key_function=None):
        super(Fasta, self).__init__(filename)
        self.key_function = key_function

        if index:
            self.index = SeqIO.index(filename, "fasta",
                    key_function=key_function)
        else:
            # SeqIO.to_dict expects a different key_function that operates on
            # the SeqRecord instead of the raw string
            _key_function = (lambda rec: key_function(rec.description)) if \
                    key_function else None 
            self.index = SeqIO.to_dict(SeqIO.parse(open(filename), "fasta"),
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
        if key in self.index:
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

    def itersizes_ordered(self):
        for rec in SeqIO.parse(open(self.filename), "fasta"):
            yield rec.name, len(rec)

    def iter_annotations(self):
        for rec in SeqIO.parse(open(self.filename), "fasta"):
            name, description = rec.name, rec.description
            yield name, description
    
    @classmethod
    def subseq(cls, fasta, start=None, stop=None, strand=None):
        """
        Take Bio.SeqRecord and slice "start:stop" from it, does proper index and
        error handling
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
        ('extract', 'given fasta file and an seq id, retrieve the sequence ' + \
                    'in fasta format'),
        ('summary', "report the real no of bases and N's in fastafiles"),
        ('uniq', 'remove records that are the same'),
        ('format', 'trim accession id to the first space'),
        ('random', 'random take some records'),
        ('diff', 'check if two FASTA records contain same information'),
        ('trim', 'given a cross_match screened fasta, trim the sequence'),
        ('pair', 'sort paired reads to .pairs.fasta and remaining to .fragments.fasta'),
        ('fastq', 'combine fasta and qual to create fastq file'),
        ('some', 'include or exclude a list of records (also performs on' + \
                 '.qual file)'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def summary(args):
    """
    %prog summary *.fasta

    Report real bases and N's in fastafiles in a tabular report
    """
    p = OptionParser(summary.__doc__)
    p.add_option("--suffix", dest="suffix", default="Mb",
            help="make the base pair couns human readable [default: %default]")
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
            help="If input reads are pairs, add trailing /1 and /2 [default: %default]")
    p.add_option("--gb", dest="gb", default=False, action="store_true",
            help="if Genbank ID, get the accession [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    infasta, outfasta = args
    gb = opts.gb
    pairs = opts.pairs

    fw = sys.stdout if outfasta=="stdout" else open(outfasta, "w")
    for i, rec in enumerate(SeqIO.parse(infasta, "fasta")):
        rec.description = ""
        if gb:
            # gi|262233616|gb|GU123895.1| Coffea arabica clone BAC
            atoms = rec.id.split("|")
            if len(atoms) >= 4:
                rec.id = atoms[3]
        if pairs:
            id = "/1" if (i % 2 == 0) else "/2"
            #rec.id = rec.id.rsplit("_", 1)[0]
            rec.id += id 
            
        SeqIO.write(rec, fw, "fasta")

    fw.close()


def get_first_rec(fastafile):
    """
    Returns the first record in the fastafile
    """
    f = list(SeqIO.parse(fastafile, "fasta"))

    if len(f) > 1:
        logging.debug("%d records found in %s, using a random one" % \
                (len(f), fastafile))

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

        if (ignore_N and 'N' in (a, b)): continue
        if a != b: break 

    if i+1==asize and i+1==bsize: 
        print_green("Two sequences match")
        match = True
    else:
        print_red("Two sequences do not match")

        snippet_size = 20 # show the context of the difference

        print_red("Sequence start to differ at position %d:" % (i+1))

        begin = max(i-snippet_size, 0)
        aend = min(i+snippet_size, asize)
        bend = min(i+snippet_size, bsize)
        
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
    p.add_option("--ignore_case", dest="ignore_case", default=False, action="store_true",
            help="ignore case when comparing sequences")
    p.add_option("--ignore_N", dest="ignore_N", default=False, action="store_true",
            help="ignore N's when comparing sequences")
    p.add_option("--rc", dest="rc", default=False, action="store_true",
            help="also consider reverse complement")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    afasta, bfasta = args
    arec = get_first_rec(afasta)
    brec = get_first_rec(bfasta)

    print tabular((arec, brec))

    asize, bsize = len(arec), len(brec)

    if asize==bsize:
        print_green("Two sequence size match (%d)" % asize)
    else:
        print_red("Two sequence size do not match (%d, %d)" % (asize, bsize))
    
    # print out the first place the two sequences diff
    print_first_difference(arec, brec, ignore_case=opts.ignore_case,
            ignore_N=opts.ignore_N, rc=opts.rc)


def some(args):
    """
    %prog some fastafile listfile outfastafile

    generate a subset of fastafile, based on a list
    """
    p = OptionParser(some.__doc__)
    p.add_option("--exclude", dest="exclude", default=False, action="store_true",
            help="output sequences not in the list file")

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(p.print_help())

    fastafile, listfile, outfastafile = args 
    outfastahandle = open(outfastafile, "w")
    qualfile = fastafile + ".qual"
    hasqual = op.exists(qualfile)
    
    names = set(x.strip() for x in open(listfile))
    if hasqual:
        outqualfile = outfastafile + ".qual"
        outqualhandle = open(outqualfile, "w")

    parser = iter_fasta_qual(fastafile, qualfile)

    num_records = 0
    for rec in parser:
        name = rec.id
        if opts.exclude:
            if name in names: continue
        else:
            if name not in names: continue
        
        rec.description = ""
        SeqIO.write([rec], outfastahandle, "fasta")
        if hasqual:
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

    fastafile = args[0]
    qualfile = fastafile + ".qual"
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

    generate .pairs.fasta and .fragments.fasta by matching records with /1/2
    into the pairs and the rest go to fragments
    """
    p = OptionParser(pair.__doc__)
    p.add_option("-d", dest="separator", 
            help="separater in the name field to reduce to the same clone " +\
                 "[e.g. GFNQ33242/1 use /, BOT01-2453H.b1 use .]")
    p.add_option("-m", dest="matepairs", default=False, action="store_true",
            help="generate .matepairs file [often used for Celera Assembler]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    fastafile = args[0]
    qualfile = fastafile + ".qual"
    hasqual = op.exists(qualfile)

    prefix = fastafile.rsplit(".", 1)[0]
    pairsfile = prefix + ".pairs.fasta"
    fragsfile = prefix + ".frags.fasta"
    pairsfw = open(pairsfile, "w")
    fragsfw = open(fragsfile, "w")

    #TODO: need a class to handle coupled fasta and qual iterating and indexing
    if opts.matepairs:
        matepairsfile = prefix + ".matepairs"
        matepairsfw = open(matepairsfile, "w")

    if hasqual:
        pairsqualfile = pairsfile + ".qual"
        pairsqualhandle = open(pairsqualfile, "w")
        fragsqualfile = fragsfile + ".qual"
        fragsqualhandle = open(fragsqualfile, "w")

    f = Fasta(args[0])
    q = SeqIO.index(qualfile, "qual")

    all_keys = list(f.iterkeys())
    all_keys.sort()
    sep = opts.separator
    for key, variants in groupby(all_keys, key=lambda x: x.split(sep, 1)[0]):
        variants = list(variants)
        paired = (len(variants) == 2)
        
        if paired and opts.matepairs:
            print >> matepairsfw, "\t".join(("%s/1" % key, "%s/2" % key))

        fw = pairsfw if paired else fragsfw
        if hasqual:
            qualfw = pairsqualhandle if paired else fragsqualhandle

        for i, var in enumerate(variants):
            rec = f[var]
            recqual = q[var]
            newid = "%s/%d" % (key, i+1)

            rec.id = newid
            rec.description = ""
            SeqIO.write([rec], fw, "fasta")
            if hasqual:
                recqual.id = newid
                recqual.description = ""
                SeqIO.write([recqual], qualfw, "qual")

    logging.debug("sequences written to `%s` and `%s`" % (pairsfile, fragsfile))
    if opts.matepairs:
        logging.debug("mates written to `%s`" % matepairsfile)


def extract(args):
    """
    %prog extract fasta query 
    
    extract query out of fasta file, query needs to be in the form of
    "seqname", or "seqname:start-stop", or "seqname:start-stop:-"
    """
    p = OptionParser(extract.__doc__)
    p.add_option('--include', dest="include", default=False, action="store_true",
            help="search full description line for match, use 'all' for all " +\
            "records [default: %default]")
    p.add_option('--exclude', dest="exclude", default=False, action="store_true",
            help="inverse search, exclude description that matches [default: %default]")

    opts, args = p.parse_args(args)

    try:
        fastafile, query = args
    except Exception, e:
        sys.exit(p.print_help())

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
            logging.error(str(e))
            sys.exit(p.print_help())

        feature["start"] = start
        feature["stop"] = stop
    else:
        start, stop = None, None

    assert start < stop or None in (start, stop), \
            "start must be less than stop, you have ({0}, {1})".format(start, stop)
    feature["strand"] = strand

    include, exclude = opts.include, opts.exclude
    # conflicting options, cannot be true at the same time
    assert not (include and exclude), "--include and --exclude cannot be "\
            "on at the same time"

    f = Fasta(fastafile)

    if include or exclude:
        for k in f.keys():
            if include and key!="all" and key not in k: continue
            if exclude and key in k: continue
            
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
            logging.error(str(e))
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
    %prog uniq fasta > uniq.fasta

    remove fasta records that are the same
    """
    p = OptionParser(uniq.__doc__)
    p.add_option("-t", "--trimname", dest="trimname",
            action="store_true", default=False,
            help="turn on the defline trim to first space [default: %default]")

    opts, args = p.parse_args(args)
    try:
        fastafile = args[0]
    except Exception as e:
        logging.error(str(e))
        sys.exit(p.print_help())

    data = {}
    for rec in _uniq_rec(fastafile):
        if opts.trimname: 
            rec.description = ""
        SeqIO.write([rec], sys.stdout, "fasta")


def random(args):
    """
    %prog random fasta 100 > random100.fasta

    take number of records randomly from fasta
    """
    p = OptionParser(random.__doc__)
    opts, args = p.parse_args(args)
    try:
        fastafile, N = args
        N = int(N)
    except Exception as e:
        logging.error(str(e))
        sys.exit(p.print_help())

    f = Fasta(fastafile)
    fw = sys.stdout
    
    for key in sample(f.keys(), N):
        rec = f[key]
        SeqIO.write([rec], fw, "fasta")


XQUAL = -1000 # default quality for X
NQUAL = 5 # default quality value for N 
QUAL = 10 # default quality value
OKQUAL = 15

def modify_qual(rec):
    qv = rec.letter_annotations['phred_quality']
    for i, (s, q) in enumerate(zip(rec.seq, qv)):
        if s=='X' or s=='x': qv[i] = XQUAL
        if s=='N' or s=='x': qv[i] = NQUAL
    return rec


def iter_fasta_qual(fastafile, qualfile, defaultqual=OKQUAL, modify=False):
    """
    used by trim, emits one SeqRecord with quality values in it
    """
    fastahandle = SeqIO.parse(fastafile, "fasta")

    if op.exists(qualfile):
        logging.debug("qual file `%s` found" % qualfile)
        qualityhandle = SeqIO.parse(qualfile, "qual")
        for rec, rec_qual in zip(fastahandle, qualityhandle):
            assert len(rec) == len(rec_qual)
            rec.letter_annotations['phred_quality'] = \
                rec_qual.letter_annotations['phred_quality']
            yield rec if not modify else modify_qual(rec)

    else:
        logging.warning("qual file `%s` not found, assume qual (%d)" % \
                (qualfile, defaultqual))
        for rec in fastahandle:
            rec.letter_annotations['phred_quality'] = [defaultqual] * len(rec) 
            yield rec if not modify else modify_qual(rec)


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
    p.add_option("-c", dest="min_length", type="int", default=15,
            help="minimum sequence length after trimming")
    p.add_option("-s", dest="score", default=QUAL,
            help="quality trimming cutoff [default: %default]")
    opts, args = p.parse_args(args)
    
    if len(args) != 2:
        sys.exit(p.print_help())

    fastafile, newfastafile = args[0:2]
    qual_suffix = ".qual"
    qualfile = fastafile + qual_suffix
    newqualfile = newfastafile + qual_suffix 

    logging.debug("Trim bad sequence from fasta file `%s` to `%s`" % \
            (fastafile, newfastafile))

    fw = open(newfastafile, "w")
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
            rec = rec[trim_start:trim_end+1]
        
        SeqIO.write([rec], fw, "fasta")
        SeqIO.write([rec], fw_qual, "qual")

    print >>sys.stderr, "A total of %d sequences modified." % trimmed
    print >>sys.stderr, "A total of %d sequences dropped (length < %d)." % \
        (dropped, opts.min_length)
            

if __name__ == '__main__':
    main()
