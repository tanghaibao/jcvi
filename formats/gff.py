#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import sys
import os
import os.path as op
import itertools
from optparse import OptionParser

from jcvi.formats.fasta import Fasta
from jcvi.formats.bed import Bed, BedLine
from jcvi.apps.base import ActionDispatcher


def main():

    actions = (
        ('bed', 'parse gff and produce bed file for particular feature type'),
        ('load', 'extract the feature (i.e. CDS) sequences and concatenate'),
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
    p.add_option("--feature", dest="feature", default="gene",
            help="the feature type to extract")

    opts, args = p.parse_args(args)
    if len(args)!=1:
        sys.exit(p.print_help())

    gff_file = args[0]
    fp = file(gff_file)

    b = Bed() 
    for row in fp:
        #chr06_pseudomolecule_IMGAG_V3   .       gene    54517   55113   .       -       .       ID=Medtr6g005000
        atoms = row.split("\t")
        chr, feature, start, stop, name = atoms[0], atoms[2], \
                atoms[3], atoms[4], atoms[-1]

        if feature!=opts.feature: continue
        
        # often the names need to have slight transformation
        chr = chr.split("_")[0]
        name = name.split(";")[0].split("=")[1].strip()

        b.append(BedLine("\t".join((chr, start, stop, name))))
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

