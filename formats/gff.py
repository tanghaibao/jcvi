#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import sys
import os
import os.path as op
import itertools
import logging

from collections import defaultdict
from urlparse import unquote
from optparse import OptionParser

from jcvi.formats.base import LineFile, must_open
from jcvi.formats.fasta import Fasta, SeqIO
from jcvi.formats.bed import Bed, BedLine
from jcvi.utils.iter import flatten
from jcvi.utils.orderedcollections import DefaultOrderedDict, parse_qs
from jcvi.apps.base import ActionDispatcher, set_outfile, mkdir, need_update, sh


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
        self.attributes_text = unquote(args[8].strip())
        self.gff3 = gff3 = "=" in self.attributes_text
        self.attributes = make_attributes(self.attributes_text, gff3=gff3)
        # key is not in the gff3 field, this indicates the conversion to accn
        self.key = key  # usually it's `ID=xxxxx;`

    def __getitem__(self, key):
        return getattr(self, key)

    def __str__(self):
        return "\t".join(str(x) for x in (self.seqid, self.source, self.type,
                self.start, self.end, self.score, self.strand, self.phase,
                self.attributes_text))

    def update_attributes(self, gff3=None):
        attributes = []
        if gff3 is None:
            gff3 = self.gff3

        sep = ";" if gff3 else "; "
        for tag, val in self.attributes.items():
            val = ",".join(val)
            val = "\"{0}\"".format(val) if " " in val or (not gff3) else val
            equal = "=" if gff3 else " "
            attributes.append(equal.join((tag, val)))

        self.attributes_text = sep.join(attributes)


    @property
    def accn(self):
        if self.key and self.key in self.attributes:
            a = self.attributes[self.key]
        else:
            a = self.attributes_text.split()
        return ",".join(a)

    id = accn

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

    def __init__(self, filename, key="ID"):
        super(Gff, self).__init__(filename)
        self.key = key

    def __iter__(self):
        fp = must_open(self.filename)
        for row in fp:
            row = row.strip()
            if row.strip() == "":
                continue
            if row[0] == '#':
                if row == FastaTag:
                    break
                continue
            yield GffLine(row, key=self.key)

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
        d = parse_qs(s)

    else:
        attributes = s.split("; ")
        d = DefaultOrderedDict(list)
        for a in attributes:
            key, val = a.strip().split(' ', 1)
            val = val.replace('"', '')
            d[key].append(val)

    for key, val in d.items():
        d[key] = list(flatten([v.split(",") for v in val]))

    return d


def main():

    actions = (
        ('bed', 'parse gff and produce bed file for particular feature type'),
        ('bed12', 'produce bed12 file for coding features'),
        ('fromgtf', 'convert gtf to gff3 format'),
        ('gtf', 'convert gff3 to gtf format'),
        ('sort', 'sort the gff file'),
        ('filter', 'filter the gff file based on Identity and Coverage'),
        ('format', 'format the gff file, change seqid, etc.'),
        ('rename', 'change the IDs within the gff3'),
        ('uniq', 'remove the redundant gene models'),
        ('liftover', 'adjust gff coordinates based on tile number'),
        ('note', 'extract certain attribute field for each feature'),
        ('load', 'extract the feature (e.g. CDS) sequences and concatenate'),
        ('extract', 'extract a particular contig from the gff file'),
        ('split', 'split the gff into one contig per file'),
        ('merge', 'merge several gff files into one'),
        ('parents', 'find the parents given a list of IDs'),
        ('children', 'find all children that belongs to the same parent'),
        ('fromgb', 'convert from gb format to gff3'),
        ('frombed', 'convert from bed format to gff3'),
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def rename(args):
    """
    %prog rename in.gff3 switch.ids > reindexed.gff3

    Change the IDs within the gff3.
    """
    from jcvi.formats.base import DictFile

    p = OptionParser(rename.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ingff3, switch = args
    switch = DictFile(switch)

    gff = Gff(ingff3)
    for g in gff:
        id, = g.attributes["ID"]
        newname = switch.get(id, id)
        g.attributes["ID"] = [newname]

        if "Parent" in g.attributes:
            parents = g.attributes["Parent"]
            g.attributes["Parent"] = [switch.get(x, x) for x in parents]

        g.update_attributes()
        print g


def parents(args):
    """
    %prog parents gffile models.ids

    Find the parents given a list of IDs in "models.ids".
    """
    p = OptionParser(parents.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gff_file, idsfile = args
    g = make_index(gff_file)
    fp = open(idsfile)
    for row in fp:
        cid = row.strip()
        b = g.parents(cid, 1).next()
        print "\t".join((cid, b.id))


def filter(args):
    """
    %prog filter gffile > filtered.gff

    Filter the gff file based on Identity and coverage. You can get this type of
    gff by using gmap:

    $ gmap -f 2
    """
    p = OptionParser(filter.__doc__)
    p.add_option("--id", default=95, type="float",
                 help="Minimum identity [default: %default]")
    p.add_option("--coverage", default=90, type="float",
                 help="Minimum coverage [default: %default]")
    p.add_option("--type", default="mRNA",
                 help="The feature to scan for the attributes [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args

    gff = Gff(gffile)
    bad = set()
    relatives = set()
    for g in gff:
        if g.type != opts.type:
            continue
        identity = float(g.attributes["Identity"][0])
        coverage = float(g.attributes["Coverage"][0])
        if identity < opts.id or coverage < opts.coverage:
            bad.add(g.accn)
            relatives.add(g.attributes["Parent"][0])

    logging.debug("{0} bad accns marked.".format(len(bad)))

    for g in gff:
        if "Parent" in g.attributes and g.attributes["Parent"][0] in bad:
            relatives.add(g.accn)

    logging.debug("{0} bad relatives marked.".format(len(relatives)))

    for g in gff:
        if g.accn in bad or g.accn in relatives:
            continue
        print g


def fix_gsac(g, notes):
    a = g.attributes

    if g.type == "gene":
        note = a["Name"]
    elif g.type == "mRNA":
        parent = a["Parent"][0]
        note = notes[parent]
    else:
        return

    a["Name"] = a["ID"]
    a["Note"] = note
    g.update_attributes()


def format(args):
    """
    %prog format gffile > formatted.gff

    Read in the gff and print it out, changing seqid, etc.
    """
    from jcvi.formats.base import DictFile

    p = OptionParser(format.__doc__)
    p.add_option("--unique", default=False, action="store_true",
                 help="Make IDs unique [default: %default]")
    p.add_option("--gff3", default=False, action="store_true",
                 help="Force to write gff3 attributes [default: %default]")
    p.add_option("--switch", help="Switch seqid from two-column file [default: %default]")
    p.add_option("--multiparents", default=False, action="store_true",
                 help="Separate features with multiple parents [default: %default]")
    p.add_option("--gsac", default=False, action="store_true",
                 help="Fix GSAC attributes [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    mapfile = opts.switch
    unique = opts.unique
    gsac = opts.gsac
    if gsac:  # setting gsac will force IDs to be unique
        unique = True

    if mapfile:
        mapping = DictFile(mapfile, delimiter="\t")

    if unique:
        dupcounts = defaultdict(int)
        gff = Gff(gffile)
        for g in gff:
            id = g.accn
            dupcounts[id] += 1
        seen = defaultdict(int)

    gff = Gff(gffile)
    notes = {}
    for g in gff:
        origid = g.seqid
        if mapfile:
            if origid in mapping:
                g.seqid = mapping[origid]
            else:
                logging.error("{0} not found in `{1}`. ID unchanged.".\
                        format(origid, mapfile))

        id = g.accn
        if unique:
            if dupcounts[id] > 1:
                seen[id] += 1
                id = "{0}-{1}".format(id, seen[id])
                g.attributes["ID"] = [id]
                g.update_attributes(gff3=True)

        if gsac and g.type == "gene":
            notes[g.accn] = g.attributes["Name"]

        pp = g.attributes.get("Parent", [])
        if opts.multiparents and len(pp) > 1:  # separate multiple parents
            id = g.attributes["ID"][0]
            for i, parent in enumerate(pp):
                g.attributes["ID"] = ["{0}-{1}".format(id, i + 1)]
                g.attributes["Parent"] = [parent]
                g.update_attributes()
                if gsac:
                    fix_gsac(g, notes)
                print g
        else:
            if opts.gff3:
                g.update_attributes(gff3=True)
            if gsac:
                fix_gsac(g, notes)
            print g


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


def get_piles(allgenes):
    """
    Before running uniq, we need to compute all the piles. The piles are a set
    of redundant features we want to get rid of. Input are a list of GffLines
    features. Output are list of list of features distinct "piles".
    """
    from jcvi.utils.range import Range, range_piles

    ranges = [Range(a.seqid, a.start, a.end, 0, i) \
                    for i, a in enumerate(allgenes)]

    for pile in range_piles(ranges):
        yield [allgenes[x] for x in pile]


def uniq(args):
    """
    %prog uniq gffile > uniq.gff

    Remove redundant gene models. For overlapping gene models, take the longest
    gene. A second scan takes only the genes selected.

    --mode controls whether you want larger feature, or higher scoring feature.
    --best controls how many redundant features to keep, e.g. 10 for est2genome.
    """
    supported_modes = ("span", "score")
    p = OptionParser(uniq.__doc__)
    p.add_option("--type", default="gene",
                 help="Types of features to non-redundify [default: %default]")
    p.add_option("--mode", default="span", choices=supported_modes,
                 help="Pile mode, one of {0} [default: %default]".\
                      format("|".join(supported_modes)))
    p.add_option("--best", default=1, type="int",
                 help="Use best N features [default: %default]")
    p.add_option("--name", default=False, action="store_true",
                 help="Non-redundify Name attribute [default: %default]")
    p.add_option("--iter", default="2", choices=("1", "2"),
                 help="Number of iterations to grab children, use 1 or 2 "\
                      "[default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    gff = Gff(gffile)
    mode = opts.mode
    bestn = opts.best
    allgenes = []
    for g in gff:
        if g.type != opts.type:
            continue
        allgenes.append(g)

    logging.debug("A total of {0} genes imported.".format(len(allgenes)))
    allgenes.sort(key=lambda x: (x.seqid, x.start))

    g = get_piles(allgenes)

    bestids = set()
    for group in g:
        if mode == "span":
            scores_group = [(- x.span, x) for x in group]
        else:
            scores_group = [(- float(x.score), x) for x in group]

        scores_group.sort()
        seen = set()
        for score, x in scores_group:
            if len(seen) >= bestn:
                break

            name = x.attributes["Name"][0] if opts.name else x.accn
            if name in seen:
                continue

            seen.add(name)
            bestids.add(x.accn)

    logging.debug("A total of {0} genes selected.".format(len(bestids)))
    logging.debug("Populate children. Iteration 1..")
    gff = Gff(gffile)
    children = set()
    for g in gff:
        if "Parent" not in g.attributes:
            continue
        for parent in g.attributes["Parent"]:
            if parent in bestids:
                children.add(g.accn)

    if opts.iter == "2":
        logging.debug("Populate children. Iteration 2..")
        gff = Gff(gffile)
        for g in gff:
            if "Parent" not in g.attributes:
                continue
            for parent in g.attributes["Parent"]:
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


def fromgtf(args):
    """
    %prog fromgtf gtffile

    Convert gtf to gff file. In gtf, the "transcript_id" will convert to "ID=",
    the "transcript_id" in exon/CDS feature will be converted to "Parent=".
    """
    p = OptionParser(fromgtf.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gtffile, = args
    gff = Gff(gtffile)
    for g in gff:
        if g.type == "transcript":
            g.type = "mRNA"
            g.attributes["ID"] = g.attributes["transcript_id"]
            g.attributes["Parent"] = g.attributes["gene_id"]
            del g.attributes["transcript_id"]
            del g.attributes["gene_id"]
        elif g.type in ("exon", "CDS"):
            g.attributes["Parent"] = g.attributes["transcript_id"]
            del g.attributes["transcript_id"]
        elif g.type == "gene":
            g.attributes["Parent"] = g.attributes["gene_id"]
        else:
            assert 0, "Doesn't know how to deal with {0}".format(g.type)

        g.update_attributes(gff3=True)
        print g


def frombed(args):
    """
    %prog frombed bed_file [--options] > gff_file

    Convert bed to gff file. In bed, the accn will convert to key='ID'
    Default type will be `match` and default source will be `source`
    """
    p = OptionParser(frombed.__doc__)
    p.add_option("--type", default="match",
                 help="GFF feature type [default: %default]")
    p.add_option("--source", default="default",
                help="GFF source qualifier [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    bed = Bed(bedfile)

    for b in bed:
        print b.gffline(type=opts.type, source=opts.source)


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
            transcript_id = g.attributes["Parent"]
        except IndexError:
            transcript_id = g.attributes["mRNA"]

        for tid in transcript_id:
            gene_id = transcript_to_gene[tid]
            g.attributes = dict(gene_id=[gene_id], transcript_id=[tid])
            g.update_attributes()

            print g


def merge(args):
    """
    %prog merge gffiles

    Merge several gff files into one. When only one file is given, it is assumed
    to be a file with a list of gff files.
    """
    p = OptionParser(merge.__doc__)
    set_outfile(p)

    opts, args = p.parse_args(args)

    nargs = len(args)
    if nargs < 1:
        sys.exit(not p.print_help())

    if nargs == 1:
        listfile, = args
        fp = open(listfile)
        gffiles = [x.strip() for x in fp]
    else:
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
    %prog extract gffile

    --contigs: Extract particular contig(s) from the gff file. If multiple contigs are
    involved, use "," to separate, e.g. "contig_12,contig_150"
    --names: Provide a file with IDs, one each line
    """
    p = OptionParser(extract.__doc__)
    p.add_option("--contigs",
                help="Extract features from certain contigs [default: %default]")
    p.add_option("--names",
                help="Extract features with certain names [default: %default]")
    p.add_option("--fasta", default=False, action="store_true",
                help="Write FASTA if available [default: %default]")
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    contigID = opts.contigs
    namesfile = opts.names

    contigID = set(contigID.split(",")) if contigID else None
    names = set(x.strip() for x in open(namesfile)) if namesfile else None

    outfile = opts.outfile
    fp = open(gffile)
    fw = must_open(outfile, "w")
    for row in fp:
        atoms = row.split()
        if len(atoms) == 0:
            continue
        tag = atoms[0]
        if row[0] == "#":
            if not (tag == RegionTag and contigID and atoms[1] not in contigID):
                print >> fw, row.rstrip()
            if tag == FastaTag:
                break
            continue

        b = GffLine(row)
        is_right_contig = (contigID and tag in contigID) or (not contigID)
        is_right_names = (names and b.attributes["Name"][0] in names) or \
                         (not names)

        if is_right_contig and is_right_names:
            print >> fw, row.rstrip()

    if not opts.fasta:
        return

    f = Fasta(gffile)
    for s in contigID:
        if s in f:
            SeqIO.write([f[s]], fw, "fasta")


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

    gff = Gff(gffile)
    seen = set()
    for g in gff:
        if attrib in g.attributes:
            keyval = (g.attributes[key][0], g.attributes[attrib][0])
            if keyval not in seen:
                print "\t".join(keyval)
                seen.add(keyval)


def bed(args):
    '''
    %prog bed gff_file [--options]

    Parses the start, stop locations of the selected features out of GFF and
    generate a bed file
    '''
    p = OptionParser(bed.__doc__)
    p.add_option("--type", dest="type", default="gene",
            help="Feature type to extract, use comma for multiple [default: %default]")
    p.add_option("--key", dest="key", default="ID",
            help="Key in the attributes to extract [default: %default]")
    set_outfile(p)

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    key = opts.key
    if key == "None":
        key = None

    type = set(x.strip() for x in opts.type.split(","))

    gff = Gff(gffile, key=key)
    b = Bed()

    for g in gff:
        if g.type not in type:
            continue

        b.append(g.bedline)

    b.sort(key=b.key)
    b.print_to_file(opts.outfile)


def make_index(gff_file):
    """
    Make a sqlite database for fast retrieval of features.
    """
    import GFFutils
    db_file = gff_file + ".db"

    if need_update(gff_file, db_file):
        if op.exists(db_file):
            os.remove(db_file)
        GFFutils.create_gffdb(gff_file, db_file)

    return GFFutils.GFFDB(db_file)


def get_parents(gff_file, parents):
    gff = Gff(gff_file)
    for g in gff:
        if g.type not in parents:
            continue
        yield g


def children(args):
    """
    %prog children gff_file

    Get the children that have the same parent.
    """
    p = OptionParser(children.__doc__)
    p.add_option("--parents", default="gene",
            help="list of features to extract, use comma to separate (e.g."
            "'gene,mRNA') [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gff_file, = args
    g = make_index(gff_file)
    parents = set(opts.parents.split(','))

    for feat in get_parents(gff_file, parents):

        cc = [c.id for c in g.children(feat.id, 1)]
        if len(cc) <= 1:
            continue

        print "\t".join(str(x) for x in \
                    (feat.id, feat.start, feat.stop, "|".join(cc)))


def load(args):
    '''
    %prog load gff_file fasta_file [--options]

    Parses the selected features out of GFF, with subfeatures concatenated.
    For example, to get the CDS sequences, do this::

    $ %prog load athaliana.gff athaliana.fa --parents mRNA --children CDS
    '''
    from jcvi.formats.fasta import Seq, SeqRecord

    p = OptionParser(load.__doc__)
    p.add_option("--parents", dest="parents", default="mRNA",
            help="list of features to extract, use comma to separate (e.g."
            "'gene,mRNA') [default: %default]")
    p.add_option("--children", dest="children", default="CDS",
            help="list of features to extract, use comma to separate (e.g."
            "'five_prime_UTR,CDS,three_prime_UTR') [default: %default]")
    p.add_option("--attribute",
            help="The attribute field to extract [default: %default]")
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    gff_file, fasta_file = args

    g = make_index(gff_file)
    f = Fasta(fasta_file, index=False)
    fw = must_open(opts.outfile, "w")

    parents = set(opts.parents.split(','))
    children_list = set(opts.children.split(','))
    attr = opts.attribute

    for feat in get_parents(gff_file, parents):

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

        description = ",".join(feat.attributes[attr]) \
                if attr and attr in feat.attributes else ""
        description = description.replace("\"", "")

        rec = SeqRecord(Seq(feat_seq), id=feat.id, description=description)
        SeqIO.write([rec], fw, "fasta")
        fw.flush()


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

    g = make_index(gffile)
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
