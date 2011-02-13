"""
parses tabular BLAST -m8 (-format 6 in BLAST+) format
"""

import sys
import logging

from itertools import groupby
from optparse import OptionParser

from jcvi.formats.base import LineFile
from jcvi.formats.coords import print_stats
from jcvi.apps.base import ActionDispatcher, debug
debug()


class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si')
 
    def __init__(self, sline):
        args = sline.split("\t")
        self.query = args[0]
        self.subject = args[1]
        self.pctid = float(args[2])
        self.hitlen = int(args[3])
        self.nmismatch = int(args[4])
        self.ngaps = int(args[5])
        self.qstart = int(args[6])
        self.qstop = int(args[7])
        self.sstart = int(args[8])
        self.sstop = int(args[9])
        self.evalue = float(args[10])
        self.score = float(args[11])
 
    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        return "\t".join(str(x) for x in \
                [getattr(self, attr) for attr in BlastLine.__slots__[:-4]])


class Blast (LineFile):
    """
    We can have a Blast class that loads entire file into memory, this is
    not very efficient for big files; when the BLAST file is generated
    by BLAST or BLAT, the file is already sorted, with same query grouped
    """
    def __init__(self, filename):
        super(Blast, self).__init__(filename)
        self.fp = open(filename)

    def iter_hits(self):
        self.fp.seek(0)
        for query, blines in groupby(self.fp, key=lambda x: BlastLine(x).query):
            blines = [BlastLine(x) for x in blines]
            blines.sort(key=lambda x: -x.score) # descending score
            yield query, blines 

    def iter_best_hit(self):
        self.fp.seek(0)
        for query, blines in groupby(self.fp, key=lambda x: BlastLine(x).query):
            blines = [BlastLine(x) for x in blines]
            best_hit = max(blines, key=lambda x: x.score)
            yield query, best_hit 
    
    @property
    def hits(self):
        """
        returns a dict with query => blastline
        """
        return dict(self.iter_hits())

    @property
    def best_hits(self):
        """
        returns a dict with query => best blasthit
        """
        return dict(self.iter_best_hit())


def get_stats(blastfile):
    
    from jcvi.utils.range import range_union

    logging.debug("report stats on `%s`" % blastfile)
    fp = open(blastfile)
    ref_ivs = []
    qry_ivs = []
    identicals = 0
    alignlen = 0

    for row in fp:
        c = BlastLine(row)
        qstart, qstop = c.qstart, c.qstop
        if qstart > qstop: qstart, qstop = qstop, qstart
        qry_ivs.append((c.query, qstart, qstop))

        sstart, sstop = c.sstart, c.sstop
        if sstart > sstop: sstart, sstop = sstop, sstart
        ref_ivs.append((c.subject, sstart, sstop))
        
        alen = sstop - sstart
        alignlen += alen
        identicals += c.pctid / 100. * alen 

    qrycovered = range_union(qry_ivs)
    refcovered = range_union(ref_ivs)
    id_pct = identicals * 100. / alignlen 

    return qrycovered, refcovered, id_pct 


def filter(args):
    """
    %prog filter test.blast > new.blast

    produce a new blast file and filter based on score
    """
    p = OptionParser(filter.__doc__)
    p.add_option("--score", dest="score", default=0., type="float",
            help="score cutoff [default: %default]")
    p.add_option("--pctid", dest="pctid", default=0., type="float",
            help="pctid cutoff [default: %default]")
    p.add_option("--hitlen", dest="hitlen", default=0., type="float",
            help="pctid cutoff [default: %default]")

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    fp = open(args[0])
    for row in fp:
        c = BlastLine(row)

        if c.score < opts.score: continue
        if c.pctid < opts.pctid: continue
        if c.hitlen < opts.hitlen: continue

        print row.rstrip()


def main():
    
    actions = (
        ('summary', 'provide summary on id% and cov%'),
        ('filter', 'filter BLAST file (based on e.g. score)'),
        ('best', 'get best BLAST hit'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def best(args):
    """
    %prog best blastfile

    print the best hit for each query in the blastfile
    """
    p = OptionParser(best.__doc__)

    opts, args = p.parse_args(args)

    if len(args)!=1:
        sys.exit(p.print_help())

    blastfile = args[0]
    b = Blast(blastfile)
    for q, bline in b.iter_best_hit():
        print bline


def summary(args):
    """
    %prog summary blastfile 
    
    provide summary on id% and cov%, for both query and reference
    """
    p = OptionParser(summary.__doc__)

    opts, args = p.parse_args(args)

    if len(args)==1:
        blastfile = args[0]
    else:
        sys.exit(p.print_help())

    qrycovered, refcovered, id_pct = get_stats(blastfile)
    print_stats(qrycovered, refcovered, id_pct)


if __name__ == '__main__':
    main()
