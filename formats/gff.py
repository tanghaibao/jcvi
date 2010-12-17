#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import sys
import os
import os.path as op
import itertools
import logging

from urlparse import parse_qs
from optparse import OptionParser

from jcvi.formats.base import LineFile
from jcvi.formats.fasta import Fasta
from jcvi.formats.bed import Bed, BedLine
from jcvi.apps.base import ActionDispatcher


Valid_strands = ('+', '-', '?', '.')
Valid_phases = ('0', '1', '2', '.')

class GffLine (object):
    """
    Specification here (http://www.sequenceontology.org/gff3.shtml)
    """
    def __init__(self, sline, key="ID"):
        args = sline.strip().split("\t")
        self.seqid = args[0]
        self.source = args[1]
        self.type = args[2]
        self.start = int(args[3])
        self.end = int(args[4])
        self.score = args[5]
        self.strand = args[6]
        assert self.strand in Valid_strands, \
                "strand must be one of %s" % Valid_strands
        self.phase = args[7]
        assert self.phase in Valid_phases, \
                "phase must be one of %s" % Valid_phases
        self.attributes = parse_qs(args[8])
        # key is not in the gff3 field, this indicates the conversion to accn
        self.key = key # usually it's `ID=xxxxx;`

    def __getitem__(self, key):
        return getattr(self, key)

    @property
    def accn(self):
        return self.attributes[self.key][0]

    @property
    def bedline(self):
        row = "\t".join((self.seqid, str(self.start-1), str(self.end), self.accn))
        return BedLine(row)


class Gff (LineFile):

    def __init__(self, filename):
        super(Gff, self).__init__(filename)
        
        fp = open(filename)
        for row in fp:
            if row[0]=='#': continue
            self.append(GffLine(row))


def main():

    actions = (
        ('bed', 'parse gff and produce bed file for particular feature type'),
        ('load', 'extract the feature (e.g. CDS) sequences and concatenate'),
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    '''
    %prog bed gff_file [--options]

    Parses the start, stop locations of the selected features out of GFF and
    generate a bed file
    '''
    p = OptionParser(bed.__doc__)
    p.add_option("--type", dest="type", default="gene",
            help="the feature type to extract")

    opts, args = p.parse_args(args)
    if len(args)!=1:
        sys.exit(p.print_help())

    fp = open(args[0])
    b = Bed() 

    seen = set()
    for row in fp: 

        if row[0]=='#': continue

        g = GffLine(row)
        if g.type!=opts.type: continue
        
        if g.seqid in seen:
            logging.error("duplicate name %s found" % g.seqid)

        b.append(g.bedline)

    b.sort(key=b.key)
    b.print_to_file()


def load(args):
    '''
    %prog load gff_file fasta_file [--options]

    Parses the selected features out of GFF, with subfeatures concatenated together.
    For example, to get the CDS sequences, do this::
        %prog athaliana.gff athaliana.fa --parents mRNA --children CDS
    '''
    p = OptionParser(load.__doc__)
    p.add_option("--parents", dest="parents", default="mRNA",
            help="list of features to extract, use comma to separate (e.g."
            "'gene,mRNA') [default: %default]")
    p.add_option("--children", dest="children", default="CDS",
            help="list of features to extract, use comma to separate (e.g."
            "'five_prime_UTR,CDS,three_prime_UTR') [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())
    
    try:
        import GFFutils
    except ImportError, e:
        logging.error(str(e))
        logging.error("You must install python library `GFFutils`")

    gff_file, fa_file = args
    parents, children = opts.parents, opts.children

    db_file = gff_file + ".db"

    if not op.exists(db_file):
        GFFutils.create_gffdb(gff_file, db_file)

    f = Fasta(fasta_file, index=False)
    g = GFFutils.GFFDB(db_file)

    parents = set(parents.split(','))
    parents_iter = [g.features_of_type(x) for x in parents]
    parents_list = itertools.chain(*parents_iter)
    children_list = set(children.split(','))

    for feat in parents_list:

        children = []
        for c in g.children(feat.id, 1):

            if c.featuretype not in children_list: continue
            child = f.sequence(dict(chr=c.chrom, start=c.start, stop=c.stop,
                strand=c.strand))
            children.append((child, c))

        if not children: 
            print >>sys.stderr, "[warning] %s has no children with type %s" \
                                    % (feat.id, ','.join(children_list))
            continue
        # sort children in incremental position
        children.sort(key=lambda x: x[1].start)
        # reverse children if negative strand
        if feat.strand=='-': children.reverse()
        feat_seq = ''.join(x[0] for x in children)

        print ">%s" % feat.id
        print feat_seq


if __name__ == '__main__':
    main()

