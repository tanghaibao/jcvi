"""
Wrapper for biopython Fasta, add option to parse sequence headers

when run as commandline,
$ python -m jcvi.formats.fasta test.fasta id

will retrieve the sequence of one record that contains the id (not necessarily
the whole description field)
"""

import sys
import logging

from optparse import OptionParser

from Bio import SeqIO
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

        fasta = self[f['chr']].seq
        
        start = f['start'] - 1 if 'start' in f else 0
        stop = f['stop'] if 'stop' in f else len(fasta)

        assert stop <= len(fasta), \
                ("stop (%d) needs to be <= than " + \
                "length of `%s` (%d)") % (stop, name, len(fasta)) 

        sequence = fasta[start:stop] 
        
        if f.get('strand') in (-1, '-1', '-'):
            sequence = sequence.reverse_complement() 
        
        if asstring:
            return str(sequence)
        
        return sequence


def main():
    
    actions = (
        ('extract', 'given fasta file and an seq id, retrieve the sequence ' + \
                    'in fasta format'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def extract(args):
    """
    %prog extract fasta query 
    
    extract query out of fasta file
    """
    p = OptionParser(extract.__doc__)

    opts, args = p.parse_args(args)

    try:
        fastafile, query = args
    except Exception, e:
        logging.error(str(e))
        sys.exit(p.print_help())

    f = Fasta(fastafile)
    for key in f.keys():
        if query not in key: continue
        
        rec = f[key]
        print ">%s\n%s" % (rec.description, rec.seq)


if __name__ == '__main__':
    main()
