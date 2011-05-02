"""
parses tabular BLAST -m8 (-format 6 in BLAST+) format
"""

import sys
import logging

from itertools import groupby
from collections import defaultdict
from optparse import OptionParser

import numpy as np

from jcvi.formats.base import LineFile
from jcvi.formats.coords import print_stats
from jcvi.utils.range import range_distance
from jcvi.graphics.histogram import histogram
from jcvi.apps.base import ActionDispatcher, debug
debug()


class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si', 'orientation')
 
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

        if self.sstart > self.sstop:
            self.sstart, self.sstop = self.sstop, self.sstart
            self.orientation = '-'
        else:
            self.orientation = '+'

 
    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        return "\t".join(str(x) for x in \
                [getattr(self, attr) for attr in BlastLine.__slots__[:12]])

    @property
    def bedline(self):
        return "\t".join(str(x) for x in (self.query, self.qstart-1, self.qstop))


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

    def iter_best_hit(self, N=1):
        self.fp.seek(0)
        for query, blines in groupby(self.fp, key=lambda x: BlastLine(x).query):
            blines = [BlastLine(x) for x in blines]
            blines.sort(key=lambda x: -x.score)
            for x in blines[:N]:
                yield query, x 
    
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
        if row[0]=='#': continue
        c = BlastLine(row)

        if c.score < opts.score: continue
        if c.pctid < opts.pctid: continue
        if c.hitlen < opts.hitlen: continue

        print row.rstrip()


def main():
    
    actions = (
        ('summary', 'provide summary on id% and cov%'),
        ('filter', 'filter BLAST file (based on e.g. score)'),
        ('best', 'get best BLAST hit per query'),
        ('pairs', 'print paired-end reads of BLAST tabular output'),
        ('bed', 'get bed file from blast'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed blastfile

    print out a bed file based on the coordinates in blast
    """
    p = OptionParser(bed.__doc__)
    
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    blastfile, = args
    fp = open(blastfile)
    fw = sys.stdout
    for row in fp:
        b = BlastLine(row)
        print >> fw, b.bedline 


def report_pairs(data, cutoff=0, dialect="blast", pairsfile=None,
        insertsfile=None):
    """
    This subroutine is used by the pairs function in blast.py and cas.py.
    Reports number of fragments and pairs as well as linked pairs
    """
    allowed_dialects = ("blast", "castab", "frgscf", "bed")
    BLAST, CASTAB, FRGSCF, BED = range(len(allowed_dialects))

    assert dialect in allowed_dialects
    dialect = allowed_dialects.index(dialect)

    num_fragments, num_pairs = 0, 0

    all_dist = []
    linked_dist = []
    # +- (forward-backward) is `innie`, -+ (backward-forward) is `outie`
    orientations = defaultdict(int)

    rs = lambda x: x[:-1]
    if dialect==BLAST:
        key = lambda x: rs(x.query)
    elif dialect==CASTAB:
        key = lambda x: rs(x.readname)
    elif dialect==FRGSCF:
        key = lambda x: rs(x.fragmentID)
    elif dialect==BED:
        key = lambda x: rs(x.accn)

    if pairsfile: pairsfw = open(pairsfile, "w")
    if insertsfile: insertsfw = open(insertsfile, "w")

    for pe, lines in groupby(data, key=key):   
        lines = list(lines)
        if len(lines) != 2:
            num_fragments += len(lines) 
            continue

        num_pairs += 1
        a, b = lines

        if dialect==BLAST:
            asubject, astart, astop = a.subject, a.sstart, a.sstop
            bsubject, bstart, bstop = b.subject, b.sstart, b.sstop

            aquery, bquery = a.query, b.query
            astrand, bstrand = a.orientation, b.orientation

        elif dialect==CASTAB:
            asubject, astart, astop = a.refnum, a.refstart, a.refstop
            bsubject, bstart, bstop = b.refnum, b.refstart, b.refstop
            if -1 in (astart, bstart): continue
            
            aquery, bquery = a.readname, b.readname
            astrand, bstrand = a.strand, b.strand

        elif dialect==FRGSCF:
            asubject, astart, astop = a.scaffoldID, a.begin, a.end
            bsubject, bstart, bstop = b.scaffoldID, b.begin, b.end

            aquery, bquery = a.fragmentID, b.fragmentID
            astrand, bstrand = a.orientation, b.orientation

        elif dialect==BED:
            asubject, astart, astop = a.seqid, a.start, a.end
            bsubject, bstart, bstop = b.seqid, b.start, b.end

            aquery, bquery = a.accn, b.accn
            astrand, bstrand = a.strand, b.strand

        dist, orientation = range_distance(\
                (asubject, astart, astop, astrand), 
                (bsubject, bstart, bstop, bstrand))

        if dist >= 0:
            all_dist.append((dist, orientation, aquery, bquery)) 

    # try to infer cutoff as twice the median until convergence
    if cutoff <= 0:
        dists = np.array([x[0] for x in all_dist], dtype="int")
        p0 = np.median(dists) 
        cutoff = int(2 * p0) # initial estimate
        logging.debug("Insert size cutoff set to {0}, ".format(cutoff) + 
            "use '--cutoff' to override")

    for dist, orientation, aquery, bquery in all_dist:
        if dist > cutoff: continue

        linked_dist.append(dist)
        if pairsfile:
            print >> pairsfw, "{0}\t{1}\t{2}".format(aquery, bquery, dist)
        orientations[orientation] += 1
    
    print >>sys.stderr, "%d fragments, %d pairs" % (num_fragments, num_pairs)
    num_links = len(linked_dist)
    
    linked_dist = np.array(linked_dist, dtype="int")
    p0 = np.median(linked_dist)
    linked_dist = np.sort(linked_dist)
    meandist = np.mean(linked_dist)
    stdev = np.std(linked_dist)
    p1, p2 = linked_dist[int(num_links * .025)], linked_dist[int(num_links * .975)]

    print >>sys.stderr, "%d pairs (%.1f%%) are linked (cutoff=%d)" % \
            (num_links, num_links*100./num_pairs, cutoff)

    print >>sys.stderr, "mean distance between PE: {0} +/-{1}".\
            format(meandist, stdev)
    print >>sys.stderr, "median distance between PE: {0}".format(p0)
    print >>sys.stderr, "95% distance range: {0} - {1}".format(p1, p2)
    print >>sys.stderr, "\nOrientations:"

    orientation_summary = []
    for orientation, count in sorted(orientations.items()):
        o = "{0}:{1}".format(orientation, count)
        orientation_summary.append(o)
        print >>sys.stderr, o 

    if insertsfile:
        print >>insertsfw, "\n".join(str(x) for x in linked_dist)
        insertsfw.close()
        prefix = insertsfile.rsplit(".", 1)[0]
        histogram(insertsfile, vmin=0, vmax=cutoff, xlabel="Insertsize", 
                title="{0} PE ({1}; median ins {2})".format(prefix, 
                    ", ".join(orientation_summary), p0))

    return meandist, stdev, p0, p1, p2


def pairs(args):
    """
    %prog pairs blastfile 
    
    report summary of blast tabular results, how many paired ends mapped, avg
    distance between paired ends, etc. Reads have to be in the form of
    `READNAME{/1,/2}`
    """
    p = OptionParser(pairs.__doc__)
    p.add_option("--cutoff", dest="cutoff", default=0, type="int",
            help="distance to call valid links between PE [default: %default]")
    p.add_option("--pairs", dest="pairsfile", 
            default=True, action="store_true",
            help="write valid pairs to pairsfile")
    p.add_option("--inserts", dest="insertsfile", default=True, 
            help="write insert sizes to insertsfile and plot distribution " + \
            "to insertsfile.pdf")
    opts, args = p.parse_args(args)

    if len(args)!=1:
        sys.exit(p.print_help())

    cutoff = opts.cutoff
    blastfile = args[0]

    basename = blastfile.split(".")[0]
    pairsfile = ".".join((basename, "pairs")) if opts.pairsfile else None
    insertsfile = ".".join((basename, "inserts")) if opts.insertsfile else None

    fp = open(blastfile)
    data = [BlastLine(row) for row in fp]
    data.sort(key=lambda x: x.query)

    report_pairs(data, cutoff, dialect="blast", pairsfile=pairsfile,
           insertsfile=insertsfile)

    
def best(args):
    """
    %prog best blastfile

    print the best hit for each query in the blastfile
    """
    p = OptionParser(best.__doc__)

    p.add_option("-N", dest="N", default=1, type="int",
            help="get best N hits [default: %default]")
    opts, args = p.parse_args(args)

    if len(args)!=1:
        sys.exit(p.print_help())

    blastfile = args[0]
    b = Blast(blastfile)
    for q, bline in b.iter_best_hit(N=opts.N):
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
