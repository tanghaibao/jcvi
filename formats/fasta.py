"""
Wrapper for biopython Fasta, add option to parse sequence headers

when run as commandline,
$ python -m jcvi.formats.fasta test.fasta id

will retrieve the sequence of one record that contains the id (not necessarily
the whole description field)
"""

import sys
import logging

from random import sample
from optparse import OptionParser

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
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def extract(args):
    """
    %prog extract fasta query 
    
    extract query out of fasta file, query needs to be in the form of
    "seqname", or "seqname:start-stop"
    """
    p = OptionParser(extract.__doc__)
    p.add_option('--all', dest="all", default=False, action="store_true",
            help="search full description line for match [default: %default]")

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

    if opts.all:
        for k in f.keys():
            if key not in k: continue
            
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
            help="turn on the defline trim to first space"
            )

    opts, args = p.parse_args(args)
    try:
        fastafile = args[0]
    except Exception as e:
        logging.error(str(e))
        sys.exit(p.print_help())

    data = {}
    for rec in _uniq_rec(fastafile):
        if opts.trimname: rec.description = ""
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


if __name__ == '__main__':
    main()
