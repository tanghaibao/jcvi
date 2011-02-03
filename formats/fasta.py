"""
Wrapper for biopython Fasta, add option to parse sequence headers
"""

import sys
import os.path as op
import logging

from random import sample
from optparse import OptionParser
from itertools import groupby

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from jcvi.formats.base import BaseFile
from jcvi.apps.base import ActionDispatcher


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

    @property
    def annotations(self):
        for rec in SeqIO.parse(open(self.filename), "fasta"):

            name, description = rec.name, rec.description
            yield name, description
    
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
        
        if 'strand' in (-1, '-1', '-'):
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
        ('uniq', 'remove records that are the same'),
        ('random', 'random take some records'),
        ('trim', 'given a cross_match screened fasta, trim the sequence'),
        ('pair', 'sort paired reads to .pairs.fasta and remaining to .fragments.fasta'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def pair(args):
    """
    %prog pair fastafile

    generate .pairs.fasta and .fragments.fasta by matching records with /1/2
    into the pairs and the rest go to fragments
    """
    p = OptionParser(pair.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    fastafile = args[0]
    prefix = fastafile.rsplit(".", 1)[0]
    pairsfile = prefix + ".pairs.fasta"
    fragsfile = prefix + ".frags.fasta"
    pairsfw = open(pairsfile, "w")
    fragsfw = open(fragsfile, "w")

    f = Fasta(args[0])
    all_keys = list(f.iterkeys())
    all_keys.sort()
    for key, variants in groupby(all_keys, key=lambda x: x.rsplit('/', 1)[0]):
        variants = list(variants)
        fw = pairsfw if len(variants)==2 else fragsfw
        for i, var in enumerate(variants):
            rec = f[var]
            rec.id = "%s/%d" % (key, i+1)
            rec.description = ""
            SeqIO.write([rec], fw, "fasta")

    logging.debug("sequences written to `%s` and `%s`" % (pairsfile, fragsfile))


def extract(args):
    """
    %prog extract fasta query 
    
    extract query out of fasta file, query needs to be in the form of
    "seqname", or "seqname:start-stop"
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
        logging.error(str(e))
        sys.exit(p.print_help())

    atoms = query.split(":")
    key = atoms[0]

    pos = "" 
    if len(atoms) >= 2:
        pos = atoms[1]

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

    f = Fasta(fastafile)

    include, exclude = opts.include, opts.exclude
    # conflicting options, cannot be true at the same time
    assert not (include and exclude), "--include and --exclude cannot be "\
            "on at the same time"

    if include or exclude:
        for k in f.keys():
            if include and key!="all" and key not in k: continue
            if exclude and key in k: continue
            
            rec = f[k]
            seq = Fasta.subseq(rec, start, stop)
            newid = rec.id
            if start is not None:
                newid += ":%d-%d" % (start, stop)

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


XQUAL = -100 # default quality for X
NQUAL = 5 # default quality value for N 
QUAL = 10 # default quality value
OKQUAL = 15

def modify_qual(rec):
    qv = rec.letter_annotations['phred_quality']
    for i, (s, q) in enumerate(zip(rec.seq, qv)):
        if s=='X' or s=='x': qv[i] = XQUAL
        if s=='N' or s=='x': qv[i] = NQUAL
    return rec


def iter_fasta_qual(fastafile, qualfile, defaultqual=OKQUAL):
    """
    used by trim, emits one SeqRecord with quality values in it
    """
    fastahandle = SeqIO.parse(fastafile, "fasta")

    if op.exists(qualfile):
        logging.warning("qual file `%s` found" % qualfile)
        qualityhandle = SeqIO.parse(qualfile, "qual")
        for rec, rec_qual in zip(fastahandle, qualityhandle):
            assert len(rec) == len(rec_qual)
            rec.letter_annotations['phred_quality'] = \
                rec_qual.letter_annotations['phred_quality']
            yield modify_qual(rec)

    else:
        logging.warning("qual file `%s` not found, assume qual (%d)" % \
                (qualfile, defaultqual))
        for rec in fastahandle:
            rec.letter_annotations['phred_quality'] = [defaultqual] * len(rec) 
            yield modify_qual(rec)


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

    print >>sys.stderr, "Trim bad sequence from fasta file `%s` to `%s`" % \
            (fastafile, newfastafile)

    fw = open(newfastafile, "w")
    fw_qual = open(newqualfile, "w")

    dropped = trimmed = 0
    
    import string

    for rec in iter_fasta_qual(fastafile, qualfile):
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
