#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import sys
import os
import os.path as op
import itertools
import logging

from collections import defaultdict
from urlparse import parse_qs
from optparse import OptionParser

from jcvi.formats.base import LineFile, must_open
from jcvi.formats.fasta import Fasta, SeqIO
from jcvi.formats.bed import Bed, BedLine
from jcvi.apps.base import ActionDispatcher, set_outfile, mkdir, sh


Valid_strands = ('+', '-', '?', '.')
Valid_phases = ('0', '1', '2', '.')
FastaTag = "##FASTA"
RegionTag = "##sequence-region"


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
                "strand must be one of {0}".format(Valid_strands)
        self.phase = args[7]
        assert self.phase in Valid_phases, \
                "phase must be one of {0}".format(Valid_phases)
        self.attributes_text = args[8].strip()
        gff3 = "=" in self.attributes_text
        self.attributes = make_attributes(self.attributes_text, gff3=gff3)
        # key is not in the gff3 field, this indicates the conversion to accn
        self.key = key  # usually it's `ID=xxxxx;`

    def __getitem__(self, key):
        return getattr(self, key)

    def __str__(self):
        return "\t".join(str(x) for x in (self.seqid, self.source, self.type,
                self.start, self.end, self.score, self.strand, self.phase,
                self.attributes_text))

    @property
    def accn(self):
        if self.key and self.key in self.attributes:
            return self.attributes[self.key][0]
        return self.attributes_text.split()[0]

    @property
    def span(self):
        return self.end - self.start + 1

    @property
    def bedline(self):
        score = "1000" if self.score == '.' else self.score
        row = "\t".join((self.seqid, str(self.start - 1),
            str(self.end), self.accn, score, self.strand))
        return BedLine(row)


class Gff (LineFile):

    def __init__(self, filename):
        super(Gff, self).__init__(filename)

        fp = open(filename)
        for row in fp:
            row = row.strip()
            if row[0] == '#':
                if row == FastaTag:
                    break
                continue
            if row.strip() == "":
                continue
            self.append(GffLine(row))

    @property
    def seqids(self):
        return set(x.seqid for x in self)


def make_attributes(s, gff3=True):
    """
    In GFF3, the last column is typically:
    ID=cds00002;Parent=mRNA00002;

    In GFF2, the last column is typically:
    Gene 22240.t000374; Note "Carbonic anhydrase"
    """
    if gff3:
        return parse_qs(s)

    attributes = s.split("; ")
    d = defaultdict(list)
    for a in attributes:
        key, val = a.strip().split(' ', 1)
        val = val.replace('"', '')
        d[key].append(val)

    return d


def main():

    actions = (
        ('bed', 'parse gff and produce bed file for particular feature type'),
        ('bed12', 'produce bed12 file for coding features'),
        ('gtf', 'convert to gtf format'),
        ('sort', 'sort the gff file'),
        ('uniq', 'remove the redundant gene models'),
        ('liftover', 'adjust gff coordinates based on tile number'),
        ('script', 'parse gmap gff and produce script for sim4db to refine'),
        ('note', 'extract certain attribute field for each feature'),
        ('load', 'extract the feature (e.g. CDS) sequences and concatenate'),
        ('extract', 'extract a particular contig from the gff file'),
        ('split', 'split the gff into one contig per file'),
        ('merge', 'merge several gff files into one'),
        ('fromgb', 'convert from gb format to gff3'),
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def liftover(args):
    """
    %prog liftover gffile > liftover.gff

    Adjust gff coordinates based on tile number. For example,
    "gannotation.asmbl.000095.7" is the 8-th tile on asmbl.000095.
    """
    p = OptionParser(liftover.__doc__)
    p.add_option("--tilesize", default=50000, type="int",
                 help="The size for each tile [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    gff = Gff(gffile)
    for g in gff:
        seqid = g.seqid
        seqid, tilenum = seqid.rsplit(".", 1)
        tilenum = int(tilenum)
        g.seqid = seqid
        offset = tilenum * opts.tilesize
        g.start += offset
        g.end += offset
        print g


def uniq(args):
    """
    %prog uniq gffile > uniq.gff

    Remove redundant gene models. For overlapping gene models, take the longest
    gene. A second scan takes only the genes selected.
    """
    from jcvi.utils.iter import pairwise
    from jcvi.utils.grouper import Grouper
    from jcvi.utils.range import range_overlap

    p = OptionParser(uniq.__doc__)
    p.add_option("--type", default="gene",
                 help="Types of features to non-redundify [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    gff = Gff(gffile)
    allgenes = []
    for g in gff:
        if g.type != opts.type:
            continue
        allgenes.append(g)

    logging.debug("A total of {0} genes imported.".format(len(allgenes)))
    allgenes.sort(key=lambda x: x.start)

    g = Grouper()
    for a, b in pairwise(allgenes):
        arange = (a.seqid, a.start, a.end)
        brange = (b.seqid, b.start, b.end)
        g.join(a)
        g.join(b)
        if range_overlap(arange, brange):
            g.join(a, b)

    bestids = set()
    for group in g:
        bestgene = max(group, key=lambda x: x.span)
        bestids.add(bestgene.accn)

    logging.debug("A total of {0} genes selected.".format(len(bestids)))
    logging.debug("Populate children. Iteration 1..")
    gff = Gff(gffile)
    children = set()
    for g in gff:
        if "Parent" not in g.attributes:
            continue
        parent = g.attributes["Parent"][0]
        if parent in bestids:
            children.add(g.accn)

    logging.debug("Populate children. Iteration 2..")
    gff = Gff(gffile)
    for g in gff:
        if "Parent" not in g.attributes:
            continue
        parent = g.attributes["Parent"][0]
        if parent in children:
            children.add(g.accn)

    logging.debug("Filter gff file..")
    gff = Gff(gffile)
    seen = set()
    for g in gff:
        accn = g.accn
        if accn in seen:
            continue
        if (g.type == opts.type and accn in bestids) or (accn in children):
            seen.add(accn)
            print g


def sort(args):
    """
    %prog sort gffile

    Sort gff file.
    """
    p = OptionParser(sort.__doc__)
    p.add_option("-i", dest="inplace", default=False, action="store_true",
                 help="Sort inplace [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    sortedgff = op.basename(gffile).rsplit(".", 1)[0] + ".sorted.gff"
    if opts.inplace:
        sortedgff = gffile

    cmd = "sort -k1,1 -k4,4n {0} -o {1}".format(gffile, sortedgff)
    sh(cmd)


def fromgb(args):
    """
    %prog fromgb gbfile

    Convert from gb format to gff3.
    """
    p = OptionParser(fromgb.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gbfile, = args
    outfile = op.basename(gbfile).rsplit(".", 1)[0] + ".bp.gff"

    cmd = "bp_genbank2gff3.pl"
    cmd += " -out stdout {0}".format(gbfile)
    sh(cmd, outfile=outfile)


def gtf(args):
    """
    %prog gtf gffile

    Convert gff to gtf file. In gtf, only exon/CDS features are important. The
    first 8 columns are the same as gff, but in the attributes field, we need to
    specify "gene_id" and "transcript_id".
    """
    p = OptionParser(gtf.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    gff = Gff(gffile)
    transcript_to_gene = {}
    for g in gff:
        if g.type == "mRNA":
            if "ID" in g.attributes and "Parent" in g.attributes:
                transcript_id = g.attributes["ID"][0]
                gene_id = g.attributes["Parent"][0]
            elif "mRNA" in g.attributes and "Gene" in g.attributes:
                transcript_id = g.attributes["mRNA"][0]
                gene_id = g.attributes["Gene"][0]
            else:
                transcript_id = g.attributes["ID"][0]
                gene_id = transcript_id
            transcript_to_gene[transcript_id] = gene_id
            continue

        if g.type not in ("CDS", "exon", "start_codon", "stop_codon"):
            continue

        try:
            transcript_id = g.attributes["Parent"] [0]
        except IndexError:
            transcript_id = g.attributes["mRNA"][0]

        transcript_id = transcript_id.split(",")  # Muliple parents
        for tid in transcript_id:
            gene_id = transcript_to_gene[tid]

            g.attributes_text = "; ".join('{0} "{1}"'.format(a, b) for a, b in \
                    zip(("gene_id", "transcript_id"), (gene_id, tid)))
            g.attributes_text += ";"
            print g


def merge(args):
    """
    %prog merge gffiles

    Merge several gff files into one.
    """
    p = OptionParser(merge.__doc__)
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    gffiles = args
    outfile = opts.outfile

    deflines = set()
    fw = must_open(outfile, "w")
    fastarecs = {}
    for gffile in gffiles:
        fp = open(gffile)
        for row in fp:
            row = row.rstrip()
            if row[0] == '#':
                if row == FastaTag:
                    break
                if row in deflines:
                    continue
                else:
                    deflines.add(row)

            print >> fw, row

        f = Fasta(gffile, lazy=True)
        for key, rec in f.iteritems_ordered():
            if key in fastarecs.keys():
                continue
            fastarecs[key] = rec

    print >> fw, FastaTag
    SeqIO.write(fastarecs.values(), fw, "fasta")


def extract(args):
    """
    %prog extract gffile contigID

    Extract particular contig(s) from the gff file. If multiple contigs are
    involved, use "," to separate, e.g. "contig_12,contig_150"
    """
    p = OptionParser(extract.__doc__)
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gffile, contigID = args
    contigID = set(contigID.split(","))

    outfile = opts.outfile
    fp = open(gffile)
    fw = must_open(outfile, "w")
    for row in fp:
        atoms = row.split()
        if len(atoms) == 0:
            continue
        tag = atoms[0]
        if row[0] == "#":
            if not (tag == RegionTag and atoms[1] not in contigID):
                print >> fw, row.rstrip()
            if tag == FastaTag:
                break
        if tag in contigID:
            print >> fw, row.rstrip()

    f = Fasta(gffile)
    for s in contigID:
        if s in f:
            SeqIO.write([f[s]], fw, "fasta")

    logging.debug("Write {0} to `{1}`.".format(",".join(contigID), outfile))


def split(args):
    """
    %prog split gffile outdir

    Split the gff into one contig per file. Will also take sequences if the file
    contains FASTA sequences.
    """
    p = OptionParser(split.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gffile, outdir = args
    mkdir(outdir)

    g = Gff(gffile)
    seqids = g.seqids

    for s in seqids:
        outfile = op.join(outdir, s + ".gff")
        extract([gffile, s, "--outfile=" + outfile])


def note(args):
    """
    %prog note gffile > tabfile

    Extract certain attribute field for each feature.
    """
    p = OptionParser(note.__doc__)
    p.add_option("--key", default="Parent",
            help="The key field to extract [default: %default]")
    p.add_option("--attribute", default="Note",
            help="The attribute field to extract [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    key = opts.key
    attrib = opts.attribute

    fp = open(gffile)
    seen = set()
    for row in fp:
        if row[0] == '#':
            continue

        g = GffLine(row)
        if attrib in g.attributes:
            keyval = (g.attributes[key][0], g.attributes[attrib][0])
            if keyval not in seen:
                print "\t".join(keyval)
                seen.add(keyval)


def script(args):
    """
    %prog script gffile cdna.fasta genome.fasta

    Parse gmap gff and produce script for sim4db to refine.
    """
    p = OptionParser(script.__doc__)

    opts, args = p.parse_args(args)
    if len(args) != 3:
        sys.exit(p.print_help())

    gffile, cdnafasta, genomefasta = args
    scriptfile = gffile + ".script"
    fp = open(gffile)
    fw = open(scriptfile, "w")
    cdnas = Fasta(cdnafasta, lazy=True)
    cdnas = dict((x, i) for (i, x) in enumerate(cdnas.iterkeys_ordered()))
    genomes = Fasta(genomefasta, lazy=True)
    genomes = dict((x, i) for (i, x) in enumerate(genomes.iterkeys_ordered()))
    extra = 50000  # 50-kb region surrounding the locus
    for row in fp:
        if row[0] == '#':
            continue

        g = GffLine(row)
        if g.type != "mRNA":
            continue

        cdna = g.attributes["Name"][0]
        genome = g.seqid
        ci = cdnas[cdna]
        gi = genomes[genome]

        strand = "-r" if g.strand == "-" else "-f"
        start, end = g.start, g.end
        start = max(0, start - extra)
        end += extra
        print >> fw, "{0} -e {1} -D {2} {3} {4}"\
                .format(strand, ci, gi, start, end)


def bed(args):
    '''
    %prog bed gff_file [--options]

    Parses the start, stop locations of the selected features out of GFF and
    generate a bed file
    '''
    p = OptionParser(bed.__doc__)
    p.add_option("--type", dest="type", default="gene",
            help="the feature type to extract [default: %default]")
    p.add_option("--key", dest="key", default="ID",
            help="the key in the attributes to extract [default: %default]")

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    key = opts.key
    if key == "None":
        key = None

    fp = open(args[0])
    b = Bed()

    seen = set()
    for row in fp:

        if row[0] == '#':
            continue

        g = GffLine(row, key=key)
        if g.type != opts.type:
            continue

        if g.seqid in seen:
            logging.error("duplicate name %s found" % g.seqid)

        b.append(g.bedline)

    b.sort(key=b.key)
    b.print_to_file()


def load(args):
    '''
    %prog load gff_file fasta_file [--options]

    Parses the selected features out of GFF, with subfeatures concatenated.
    For example, to get the CDS sequences, do this::
        %prog athaliana.gff athaliana.fa --parents mRNA --children CDS
    '''
    import GFFutils

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

    gff_file, fasta_file = args
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

            if c.featuretype not in children_list:
                continue
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
        if feat.strand == '-':
            children.reverse()
        feat_seq = ''.join(x[0] for x in children)

        print ">%s" % feat.id
        print feat_seq


def bed12(args):
    """
    %prog bed12 gffile > bedfile

    Produce bed12 file for coding features. The exons will be converted to blocks.
    The CDS range will be shown between thickStart to thickEnd. For reference,
    bed format consists of the following fields:

    1. chrom
    2. chromStart
    3. chromEnd
    4. name
    5. score
    6. strand
    7. thickStart
    8. thickEnd
    9. itemRgb
    10. blockCount
    11. blockSizes
    12. blockStarts
    """
    import GFFutils

    p = OptionParser(bed12.__doc__)
    p.add_option("--parent", default="mRNA",
            help="Top feature type [default: %default]")
    p.add_option("--block", default="exon",
            help="Feature type for regular blocks [default: %default]")
    p.add_option("--thick", default="CDS",
            help="Feature type for thick blocks [default: %default]")
    set_outfile(p)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    parent, block, thick = opts.parent, opts.block, opts.thick
    outfile = opts.outfile

    dbfile = gffile + ".db"

    if not op.exists(dbfile):
        GFFutils.create_gffdb(gffile, dbfile)

    g = GFFutils.GFFDB(dbfile)
    fw = must_open(outfile, "w")

    for f in g.features_of_type(parent):

        chrom = f.chrom
        chromStart = f.start - 1
        chromEnd = f.stop
        name = f.id
        score = 0
        strand = f.strand
        thickStart = 1e15
        thickEnd = 0
        blocks = []

        for c in g.children(name, 1):

            cstart, cend = c.start - 1, c.stop

            if c.featuretype == block:
                blockStart = cstart - chromStart
                blockSize = cend - cstart
                blocks.append((blockStart, blockSize))

            elif c.featuretype == thick:
                thickStart = min(thickStart, cstart)
                thickEnd = max(thickEnd, cend)

        blocks.sort()
        blockStarts, blockSizes = zip(*blocks)
        blockCount = len(blocks)
        blockSizes = ",".join(str(x) for x in blockSizes) + ","
        blockStarts = ",".join(str(x) for x in blockStarts) + ","
        itemRgb = 0

        print >> fw, "\t".join(str(x) for x in (chrom, chromStart, chromEnd, \
                name, score, strand, thickStart, thickEnd, itemRgb,
                blockCount, blockSizes, blockStarts))


if __name__ == '__main__':
    main()
