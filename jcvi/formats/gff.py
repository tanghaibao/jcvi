#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import os
import os.path as op
import logging
import re

from collections import defaultdict
from more_itertools import flatten
from urllib.parse import quote, unquote

from jcvi.utils.cbook import AutoVivification
from jcvi.formats.base import DictFile, LineFile, must_open, is_number
from jcvi.formats.fasta import Fasta, SeqIO
from jcvi.formats.bed import Bed, BedLine, natsorted
from jcvi.annotation.reformat import atg_name
from jcvi.utils.range import range_minmax
from jcvi.utils.orderedcollections import DefaultOrderedDict, OrderedDict, parse_qs
from jcvi.apps.base import (
    OptionParser,
    OptionGroup,
    ActionDispatcher,
    mkdir,
    parse_multi_values,
    need_update,
    sh,
)


Valid_strands = ("+", "-", "?", ".")
Valid_phases = ("0", "1", "2", ".")
FastaTag = "##FASTA"
RegionTag = "##sequence-region"
valid_gff_parent_child = {
    "match": "match_part",
    "cDNA_match": "match_part",
    "EST_match": "match_part",
    "nucleotide_to_protein_match": "match_part",
    "expressed_sequence_match": "match_part",
    "protein_match": "match_part",
    "transposable_element": "transposon_fragment",
    "gene": "mRNA",
    "mRNA": "exon,CDS,five_prime_UTR,three_prime_UTR",
}
valid_gff_to_gtf_type = {
    "exon": "exon",
    "pseudogenic_exon": "exon",
    "CDS": "CDS",
    "start_codon": "start_codon",
    "stop_codon": "stop_codon",
    "five_prime_UTR": "5UTR",
    "three_prime_UTR": "3UTR",
}
valid_gff_type = tuple(valid_gff_parent_child.keys())
reserved_gff_attributes = (
    "ID",
    "Name",
    "Alias",
    "Parent",
    "Target",
    "Gap",
    "Derives_from",
    "Note",
    "Dbxref",
    "Ontology_term",
    "Is_circular",
)
multiple_gff_attributes = ("Parent", "Alias", "Dbxref", "Ontology_term")
safechars = " /:?~#+!$'@()*[]|"


class GffLine(object):
    """
    Specification here (http://www.sequenceontology.org/gff3.shtml)
    """

    def __init__(
        self,
        sline,
        key="ID",
        gff3=True,
        line_index=None,
        strict=True,
        append_source=False,
        append_ftype=False,
        score_attrib=False,
        keep_attr_order=True,
        compute_signature=False,
    ):
        sline = sline.strip()
        args = sline.split("\t")
        if len(args) != 9:
            args = sline.split()
        if strict:
            assert len(args) == 9, "Malformed line ({0} columns != 9): {1}".format(
                len(args), args
            )
        self.seqid = args[0]
        self.source = args[1]
        self.type = args[2]
        self.start = int(args[3])
        self.end = int(args[4])
        self.score = args[5]
        self.strand = args[6]
        assert self.strand in Valid_strands, "strand must be one of {0}".format(
            Valid_strands
        )
        self.phase = args[7]
        assert self.phase in Valid_phases, "phase must be one of {0}".format(
            Valid_phases
        )
        self.attributes_text = "" if len(args) <= 8 else args[8].strip()
        self.attributes = make_attributes(
            self.attributes_text, gff3=gff3, keep_attr_order=keep_attr_order
        )
        # key is not in the gff3 field, this indicates the conversion to accn
        self.key = key  # usually it's `ID=xxxxx;`
        self.gff3 = gff3

        if append_ftype and self.key in self.attributes:
            # if `append_ftype` is True, append the gff `self.type`
            # to `self.key`. use this option to enhance the `self.accn`
            # column in bed file
            self.attributes[self.key][0] = ":".join(
                (self.type, self.attributes[self.key][0])
            )

        if append_source and self.key in self.attributes:
            # if `append_source` is True, append the gff `self.source`
            # to `self.key`. use this option to enhance the `self.accn`
            # column in bed file
            self.attributes[self.key][0] = ":".join(
                (self.source, self.attributes[self.key][0])
            )

        if (
            score_attrib
            and score_attrib in self.attributes
            and is_number(self.attributes[score_attrib][0])
        ):
            # if `score_attrib` is specified, check if it is indeed an
            # attribute or not. If yes, check if the value of attribute
            # is numeric or not. If not, keep original GFF score value
            self.score = self.attributes[score_attrib][0]

        if line_index is not None and is_number(line_index):
            # if `line_index` in provided, initialize an idx variable
            # used to autcompute the ID for a feature
            self.idx = line_index

        if compute_signature:
            # if `compute_signature` is specified, compute a signature for
            # the gff line and store in variable `sign`
            self.sign = self.signature

    def __getitem__(self, key):
        return getattr(self, key)

    def __str__(self):
        return "\t".join(
            str(x)
            for x in (
                self.seqid,
                self.source,
                self.type,
                self.start,
                self.end,
                self.score,
                self.strand,
                self.phase,
                self.attributes_text,
            )
        )

    def get_attr(self, key, first=True):
        if key in self.attributes:
            if first:
                return self.attributes[key][0]
            return self.attributes[key]
        return None

    def set_attr(
        self, key, value, update=False, append=False, dbtag=None, urlquote=False
    ):
        if value is None:
            self.attributes.pop(key, None)
        else:
            if key == "Dbxref" and dbtag:
                value = value.split(",")
                value = ["{0}:{1}".format(dbtag, x) for x in value]
            if type(value) is not list:
                value = [value]
            if key not in self.attributes or not append:
                self.attributes[key] = []
            self.attributes[key].extend(value)
        if update:
            self.update_attributes(gff3=self.gff3, urlquote=urlquote)

    def update_attributes(self, skipEmpty=True, gff3=True, gtf=None, urlquote=True):
        attributes = []
        if gtf:
            gff3 = None
        elif gff3 is None:
            gff3 = self.gff3

        sep = ";" if gff3 else "; "
        for tag, val in self.attributes.items():
            if not val and skipEmpty:
                continue
            val = ",".join(val)
            val = '"{0}"'.format(val) if (" " in val and (not gff3)) or gtf else val
            equal = "=" if gff3 else " "
            if urlquote:
                sc = safechars
                if tag in multiple_gff_attributes:
                    sc += ","
                val = quote(val, safe=sc)
            attributes.append(equal.join((tag, val)))

        self.attributes_text = sep.join(attributes)
        if gtf:
            self.attributes_text += ";"

    def update_tag(self, old_tag, new_tag):
        if old_tag not in self.attributes:
            return
        self.attributes[new_tag] = self.attributes[old_tag]
        del self.attributes[old_tag]

    @property
    def accn(self):
        if self.key:  # GFF3 format
            if self.key not in self.attributes:
                a = ["{0}_{1}".format(str(self.type).lower(), self.idx)]
            else:
                a = self.attributes[self.key]
        else:  # GFF2 format
            a = self.attributes_text.split()
        return quote(",".join(a), safe=safechars)

    id = accn

    @property
    def name(self):
        return self.attributes["Name"][0] if "Name" in self.attributes else None

    @property
    def parent(self):
        return self.attributes["Parent"][0] if "Parent" in self.attributes else None

    @property
    def span(self):
        return self.end - self.start + 1

    @property
    def bedline(self):
        score = "0" if self.score == "." else self.score
        row = "\t".join(
            (
                self.seqid,
                str(self.start - 1),
                str(self.end),
                self.accn,
                score,
                self.strand,
            )
        )
        return BedLine(row)

    @property
    def signature(self):
        """
        create a unique signature for any GFF line based on joining
        columns 1,2,3,4,5,7,8 (into a comma separated string)
        """
        sig_elems = [
            self.seqid,
            self.source,
            self.type,
            self.start,
            self.end,
            self.strand,
            self.phase,
        ]
        if re.search("exon|CDS|UTR", self.type):
            parent = self.get_attr("Parent")
            if parent:
                (locus, iso) = atg_name(parent, retval="locus,iso", trimpad0=False)
                if locus:
                    sig_elems.append(locus)
        else:
            sig_elems.extend([self.accn])

        return ",".join(str(elem) for elem in sig_elems)


class Gff(LineFile):
    def __init__(
        self,
        filename,
        key="ID",
        strict=True,
        append_source=False,
        append_ftype=False,
        score_attrib=False,
        keep_attr_order=True,
        make_gff_store=False,
        compute_signature=False,
    ):
        super(Gff, self).__init__(filename)
        self.make_gff_store = make_gff_store
        self.gff3 = True
        if self.make_gff_store:
            self.gffstore = []
            gff = Gff(
                self.filename,
                key=key,
                strict=True,
                append_source=append_source,
                append_ftype=append_ftype,
                score_attrib=score_attrib,
                keep_attr_order=keep_attr_order,
                compute_signature=compute_signature,
            )
            for g in gff:
                self.gffstore.append(g)
        else:
            self.key = key
            self.strict = strict
            self.append_source = append_source
            self.append_ftype = append_ftype
            self.score_attrib = score_attrib
            self.keep_attr_order = keep_attr_order
            self.compute_signature = compute_signature
            if filename in ("-", "stdin") or filename.endswith(".gz"):
                if ".gtf" in filename:
                    self.gff3 = False
                    logging.debug("File is not gff3 standard.")
                return

            self.set_gff_type()

    def set_gff_type(self):
        # Determine file type
        row = None
        for row in self:
            break
        gff3 = False if not row else "=" in row.attributes_text
        if not gff3:
            logging.debug("File is not gff3 standard.")

        self.gff3 = gff3
        self.fp.seek(0)

    def __iter__(self):
        if self.make_gff_store:
            for row in self.gffstore:
                yield row
        else:
            self.fp = must_open(self.filename)
            for idx, row in enumerate(self.fp):
                row = row.strip()
                if row.strip() == "":
                    continue
                if row[0] == "#":
                    if row == FastaTag:
                        break
                    continue
                yield GffLine(
                    row,
                    key=self.key,
                    line_index=idx,
                    strict=self.strict,
                    append_source=self.append_source,
                    append_ftype=self.append_ftype,
                    score_attrib=self.score_attrib,
                    keep_attr_order=self.keep_attr_order,
                    compute_signature=self.compute_signature,
                    gff3=self.gff3,
                )

    @property
    def seqids(self):
        return set(x.seqid for x in self)


class GffFeatureTracker(object):
    def __init__(self):
        self.ftype = "exon|CDS|UTR|fragment"
        self.tracker = {}
        self.symbolstore = {}

    def track(self, parent, g):
        if re.search(self.ftype, g.type):
            if parent not in self.tracker:
                self.tracker[parent] = {}
            if g.type not in self.tracker[parent]:
                self.tracker[parent][g.type] = set()
            self.tracker[parent][g.type].add((g.start, g.end, g.sign))

    def _sort(self, parent, ftype, reverse=False):
        if not isinstance(self.tracker[parent][ftype], list):
            self.tracker[parent][ftype] = sorted(
                list(self.tracker[parent][ftype]),
                key=lambda x: (x[0], x[1]),
                reverse=reverse,
            )

    def feat_index(self, parent, ftype, strand, feat_tuple):
        reverse = True if strand == "-" else False
        self._sort(parent, ftype, reverse=reverse)
        return self.tracker[parent][ftype].index(feat_tuple)

    def store_symbol(self, g):
        for symbol_attr in ("symbol", "ID"):
            if symbol_attr in g.attributes:
                break
        self.symbolstore[g.accn] = g.get_attr(symbol_attr)

    def get_symbol(self, parent):
        return self.symbolstore[parent]


def make_attributes(s, gff3=True, keep_attr_order=True):
    """
    In GFF3, the last column is typically:
    ID=cds00002;Parent=mRNA00002;

    In GFF2, the last column is typically:
    Gene 22240.t000374; Note "Carbonic anhydrase"
    """
    if gff3:
        """
        hack: temporarily replace the '+' sign in the attributes column
        with the string 'PlusSign' to prevent urlparse.parse_qsl() from
        replacing the '+' sign with a space
        """
        s = s.replace("+", "PlusSign")
        d = parse_qs(s, separator=";", keep_attr_order=keep_attr_order)
        for key in d:
            d[key][0] = unquote(d[key][0].replace("PlusSign", "+").replace('"', ""))
    else:
        attributes = s.split(";")
        d = DefaultOrderedDict(list) if keep_attr_order else defaultdict(list)
        for a in attributes:
            a = a.strip()
            if " " not in a:
                continue
            key, val = a.split(" ", 1)
            val = unquote(val.replace('"', "").replace("=", " ").strip())
            d[key].append(val)

    for key, val in d.items():
        d[key] = list(flatten([v.split(",") for v in val]))

    return d


def to_range(obj, score=None, id=None, strand=None):
    """
    Given a gffutils object, convert it to a range object
    """
    from jcvi.utils.range import Range

    if score or id:
        _score = score if score else obj.score
        _id = id if id else obj.id
        return Range(
            seqid=obj.seqid, start=obj.start, end=obj.end, score=_score, id=_id
        )
    elif strand:
        return obj.seqid, obj.start, obj.end, obj.strand

    return obj.seqid, obj.start, obj.end


def main():

    actions = (
        ("addparent", "merge sister features and infer their parent"),
        ("bed", "parse gff and produce bed file for particular feature type"),
        ("bed12", "produce bed12 file for coding features"),
        ("fromgtf", "convert gtf to gff3 format"),
        ("gtf", "convert gff3 to gtf format"),
        ("gb", "convert gff3 to genbank format"),
        ("sort", "sort the gff file"),
        ("filter", "filter the gff file based on Identity and Coverage"),
        ("sizes", "calculate sizes of features in gff file"),
        ("format", "format the gff file, change seqid, etc."),
        (
            "fixboundaries",
            "fix boundaries of parent features by range chaining child features",
        ),
        ("chain", "fill in parent features by chaining children"),
        ("rename", "change the IDs within the gff3"),
        ("uniq", "remove the redundant gene models"),
        ("liftover", "adjust gff coordinates based on tile number"),
        ("note", "extract certain attribute field for each feature"),
        ("load", "extract the feature (e.g. CDS) sequences and concatenate"),
        ("extract", "extract contig or features from gff file"),
        ("split", "split the gff into one contig per file"),
        ("merge", "merge several gff files into one"),
        ("parents", "find the parents given a list of IDs"),
        ("children", "find all children that belongs to the same parent"),
        ("frombed", "convert from bed format to gff3"),
        ("fromsoap", "convert from soap format to gff3"),
        ("gapsplit", "split alignment GFF3 at gaps based on CIGAR string"),
        ("orient", "orient the coding features based on translation"),
        ("splicecov", "tag gff introns with coverage info from junctions.bed"),
        ("summary", "print summary stats for features of different types"),
        ("cluster", "cluster transcripts based on shared splicing structure"),
        (
            "fixpartials",
            "fix 5/3 prime partial transcripts, locate nearest in-frame start/stop",
        ),
    )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def addparent(args):
    """
    %prog addparent file.gff

    Merge sister features and infer parents.
    """
    p = OptionParser(addparent.__doc__)
    p.add_option("--childfeat", default="CDS", help="Type of children feature")
    p.add_option("--parentfeat", default="mRNA", help="Type of merged feature")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gff_file,) = args
    gff = Gff(gff_file)
    data = defaultdict(list)
    for g in gff:
        if g.type != opts.childfeat:
            continue
        data[g.parent].append(g)

    logging.debug(
        "A total of {0} {1} features clustered".format(len(data), opts.childfeat)
    )

    parents = []
    for parent, dd in data.items():
        d = dd[0]
        start, end = min(x.start for x in dd), max(x.end for x in dd)
        gffline = "\t".join(
            str(x)
            for x in (
                d.seqid,
                d.source,
                opts.parentfeat,
                start,
                end,
                ".",
                d.strand,
                ".",
                "ID={0};Name={0}".format(parent),
            )
        )
        parents.append(GffLine(gffline))
    parents.sort(key=lambda x: (x.seqid, x.start))
    logging.debug("Merged feature sorted")

    fw = must_open(opts.outfile, "w")
    for parent in parents:
        print(parent, file=fw)
        parent_id = parent.id
        for d in data[parent_id]:
            if d.accn == parent_id:
                new_id = "{0}.{1}1".format(parent_id, opts.childfeat)
                d.set_attr("ID", new_id)
                d.set_attr("Name", new_id, update=True)
            print(d, file=fw)
    fw.close()


def _fasta_slice(fasta, seqid, start, stop, strand):
    """
    Return slice of fasta, given (seqid, start, stop, strand)
    """
    _strand = 1 if strand == "+" else -1
    return fasta.sequence(
        {"chr": seqid, "start": start, "stop": stop, "strand": _strand}
    )


def is_valid_codon(codon, type="start"):
    """
    Given a codon sequence, check if it is a valid start/stop codon
    """
    if len(codon) != 3:
        return False

    if type == "start":
        if codon != "ATG":
            return False
    elif type == "stop":
        if not any(_codon == codon for _codon in ("TGA", "TAG", "TAA")):
            return False
    else:
        logging.error(
            "`{0}` is not a valid codon type. ".format(type)
            + "Should be one of (`start` or `stop`)"
        )
        sys.exit()

    return True


def scan_for_valid_codon(codon_span, strand, seqid, genome, type="start"):
    """
    Given a codon span, strand and reference seqid, scan upstream/downstream
    to find a valid in-frame start/stop codon
    """
    s, e = codon_span[0], codon_span[1]
    while True:
        if (type == "start" and strand == "+") or (type == "stop" and strand == "-"):
            s, e = s - 3, e - 3
        else:
            s, e = s + 3, e + 3

        codon = _fasta_slice(genome, seqid, s, e, strand)
        is_valid = is_valid_codon(codon, type=type)
        if not is_valid:
            if type == "start":
                # if we are scanning upstream for a valid start codon,
                # stop scanning when we encounter a stop
                if is_valid_codon(codon, type="stop"):
                    return None, None
            elif type == "stop":
                # if we are scanning downstream for a valid stop codon,
                # stop scanning when we encounter a start
                if is_valid_codon(codon, type="start"):
                    return None, None
            continue
        break

    return s, e


def fixpartials(args):
    """
    %prog fixpartials genes.gff genome.fasta partials.ids

    Given a gff file of features, fix partial (5'/3' incomplete) transcripts
    by trying to locate nearest in-frame start/stop codon
    """
    p = OptionParser(fixpartials.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    (
        gffile,
        gfasta,
        partials,
    ) = args

    gff = make_index(gffile)
    genome = Fasta(gfasta, index=True)
    partials = LineFile(partials, load=True).lines

    # all_transcripts = [f.id for f in gff.features_of_type("mRNA", \
    #     order_by=("seqid", "start"))]
    seen = set()
    fw = must_open(opts.outfile, "w")
    for gene in gff.features_of_type("gene", order_by=("seqid", "start")):
        children = AutoVivification()
        cflag = False
        transcripts = list(gff.children(gene, level=1, order_by="start"))
        for transcript in transcripts:
            trid, seqid, strand = transcript.id, transcript.seqid, transcript.strand

            for child in gff.children(transcript, order_by="start"):
                ftype = child.featuretype
                if ftype not in children[trid]:
                    children[trid][ftype] = []
                children[trid][ftype].append(child)

            five_prime, three_prime = True, True
            nstart, nstop = (None, None), (None, None)
            cds_span = [children[trid]["CDS"][0].start, children[trid]["CDS"][-1].stop]
            new_cds_span = [x for x in cds_span]

            start_codon = (cds_span[0], cds_span[0] + 2)
            stop_codon = (cds_span[1] - 2, cds_span[1])
            if strand == "-":
                start_codon, stop_codon = stop_codon, start_codon

            if trid in partials:
                seen.add(trid)
                start_codon_fasta = _fasta_slice(
                    genome, seqid, start_codon[0], start_codon[1], strand
                )
                stop_codon_fasta = _fasta_slice(
                    genome, seqid, stop_codon[0], stop_codon[1], strand
                )

                if not is_valid_codon(start_codon_fasta, type="start"):
                    five_prime = False
                    nstart = scan_for_valid_codon(
                        start_codon, strand, seqid, genome, type="start"
                    )

                if not is_valid_codon(stop_codon_fasta, type="stop"):
                    three_prime = False
                    nstop = scan_for_valid_codon(
                        stop_codon, strand, seqid, genome, type="stop"
                    )

                logging.debug(
                    "feature={0} ({1})".format(trid, strand)
                    + ", 5'={0}, 3'={1}".format(five_prime, three_prime)
                    + ", {0} <== {1} ==> {2}".format(
                        nstart if strand == "+" else nstop,
                        cds_span,
                        nstop if strand == "+" else nstart,
                    )
                )

            if not five_prime or not three_prime:
                if nstart != (None, None) and (start_codon != nstart):
                    i = 0 if strand == "+" else 1
                    new_cds_span[i] = nstart[i]
                if nstop != (None, None) and (stop_codon != nstop):
                    i = 1 if strand == "+" else 0
                    new_cds_span[i] = nstop[i]
                new_cds_span.sort()

            if set(cds_span) != set(new_cds_span):
                cflag = True
                # if CDS has been extended, appropriately adjust all relevent
                # child feature (CDS, exon, UTR) coordinates
                for ftype in children[trid]:
                    for idx in range(len(children[trid][ftype])):
                        child_span = (
                            children[trid][ftype][idx].start,
                            children[trid][ftype][idx].stop,
                        )
                        if ftype in ("exon", "CDS"):
                            # if exons/CDSs, adjust start and stop according to
                            # new CDS start and stop, respectively
                            if child_span[0] == cds_span[0]:
                                children[trid][ftype][idx].start = new_cds_span[0]
                            if child_span[1] == cds_span[1]:
                                children[trid][ftype][idx].stop = new_cds_span[1]
                        elif ftype.endswith("UTR"):
                            # if *_prime_UTR, adjust stop according to new CDS start and
                            #                 adjust start according to new CDS stop
                            if child_span[1] == cds_span[0]:
                                children[trid][ftype][idx].stop = new_cds_span[0]
                            if child_span[0] == cds_span[1]:
                                children[trid][ftype][idx].start = new_cds_span[1]

                transcript.start, transcript.stop = (
                    children[trid]["exon"][0].start,
                    children[trid]["exon"][-1].stop,
                )

        if cflag:
            _gene_span = range_minmax([(tr.start, tr.stop) for tr in transcripts])
            gene.start, gene.stop = _gene_span[0], _gene_span[1]

        # print gff file
        print(gene, file=fw)
        for transcript in transcripts:
            trid = transcript.id
            print(transcript, file=fw)
            for cftype in children[trid]:
                for child in children[trid][cftype]:
                    print(child, file=fw)

    fw.close()


def sizes(args):
    """
    %prog sizes gffile

    Given a gff file of features, calculate the sizes of chosen parent feature
    based on summation of sizes of child features.

    For example, for parent 'mRNA' and child 'CDS' feature types, calcuate sizes of
    mRNA by summing the sizes of the disjoint CDS parts.
    """
    p = OptionParser(sizes.__doc__)
    p.set_outfile()
    p.add_option(
        "--parents",
        dest="parents",
        default="mRNA",
        help="parent feature(s) for which size is to be calculated",
    )
    p.add_option(
        "--child",
        dest="child",
        default="CDS",
        help="child feature to use for size calculations",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    parents, cftype = set(opts.parents.split(",")), opts.child

    gff = make_index(gffile)

    fw = must_open(opts.outfile, "w")
    for parent in parents:
        for feat in gff.features_of_type(parent, order_by=("seqid", "start")):
            fsize = 0
            fsize = (
                feat.end - feat.start + 1
                if cftype == parent
                else gff.children_bp(feat, child_featuretype=cftype)
            )
            print("\t".join(str(x) for x in (feat.id, fsize)), file=fw)
        fw.close()


def cluster(args):
    """
    %prog cluster gffile

    Given a gff file of gene structures (multiple transcripts per gene locus),
    cluster/consolidate all transcripts based on shared splicing structure.

    If `slop` is enabled, clustering/consolidation will collapse any variation
    in terminal UTR lengths, keeping only the longest as representative.
    """
    from jcvi.utils.grouper import Grouper
    from itertools import combinations

    p = OptionParser(cluster.__doc__)
    p.add_option(
        "--slop",
        default=False,
        action="store_true",
        help="allow minor variation in terminal 5'/3' UTR" + " start/stop position",
    )
    p.add_option(
        "--inferUTR",
        default=False,
        action="store_true",
        help="infer presence of UTRs from exon coordinates",
    )
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    slop = opts.slop
    inferUTR = opts.inferUTR

    gff = make_index(gffile)

    fw = must_open(opts.outfile, "w")
    print("##gff-version	3", file=fw)
    seen = {}
    for gene in gff.features_of_type("gene", order_by=("seqid", "start")):
        g = Grouper()
        mrnas = list(
            combinations(
                [
                    mrna
                    for mrna in gff.children(gene, featuretype="mRNA", order_by="start")
                ],
                2,
            )
        )
        if len(mrnas) > 0:
            for mrna1, mrna2 in mrnas:
                mrna1s, mrna2s = (
                    gff.children_bp(mrna1, child_featuretype="exon"),
                    gff.children_bp(mrna2, child_featuretype="exon"),
                )
                g.join((mrna1.id, mrna1s))
                g.join((mrna2.id, mrna2s))

                if match_subfeats(mrna1, mrna2, gff, gff, featuretype="CDS"):
                    res = []
                    ftypes = (
                        ["exon"] if inferUTR else ["five_prime_UTR", "three_prime_UTR"]
                    )
                    for ftype in ftypes:
                        res.append(
                            match_subfeats(
                                mrna1, mrna2, gff, gff, featuretype=ftype, slop=slop
                            )
                        )

                    if all(res):
                        g.join((mrna1.id, mrna1s), (mrna2.id, mrna2s))
        else:
            for mrna1 in gff.children(gene, featuretype="mRNA", order_by="start"):
                mrna1s = gff.children_bp(mrna1, child_featuretype="exon")
                g.join((mrna1.id, mrna1s))

        print(gene, file=fw)
        for group in sorted(g):
            group.sort(key=lambda x: x[1], reverse=True)
            mrnas = [el[0] for el in group]
            m = mrnas[0]

            _mrnaid = []
            for x in mrnas:
                if x not in _mrnaid:
                    _mrnaid.append(x)
            mrnaid = "{0}".format("-".join(_mrnaid))
            if mrnaid not in seen:
                seen[mrnaid] = 0
            else:
                seen[mrnaid] += 1
                mrnaid = "{0}-{1}".format(mrnaid, seen[mrnaid])

            _mrna = gff[m]
            _mrna.attributes["ID"] = [mrnaid]
            _mrna.attributes["Parent"] = [gene.id]
            children = gff.children(m, order_by="start")
            print(_mrna, file=fw)
            for child in children:
                child.attributes["ID"] = ["{0}".format(child.id)]
                child.attributes["Parent"] = [mrnaid]
                print(child, file=fw)

    fw.close()


def summary(args):
    """
    %prog summary gffile

    Print summary stats for features of different types.
    """
    from jcvi.formats.base import SetFile
    from jcvi.formats.bed import BedSummary
    from jcvi.utils.table import tabulate

    p = OptionParser(summary.__doc__)
    p.add_option(
        "--isoform",
        default=False,
        action="store_true",
        help="Find longest isoform of each id",
    )
    p.add_option("--ids", help="Only include features from certain IDs")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gff_file,) = args
    ids = opts.ids

    if ids:
        ids = SetFile(ids)
        logging.debug("Total ids loaded: {0}".format(len(ids)))

        if opts.isoform:
            pids = set()
            gff = Gff(gff_file)
            for g in gff:
                if g.type != "mRNA":
                    continue
                if g.parent not in ids:
                    continue
                if "longest" not in g.attributes:
                    pids = set(x + ".1" for x in ids)
                    break
                if g.attributes["longest"][0] == "0":
                    continue
                pids.add(g.id)
            ids = pids
            logging.debug("After checking longest: {0}".format(len(ids)))

        # Collects aliases
        gff = Gff(gff_file)
        for g in gff:
            if g.name in ids:
                ids.add(g.id)
        logging.debug("Total ids including aliases: {0}".format(len(ids)))

    gff = Gff(gff_file)
    beds = defaultdict(list)
    for g in gff:
        if ids and not (g.id in ids or g.name in ids or g.parent in ids):
            continue

        beds[g.type].append(g.bedline)

    table = {}
    for type, bb in sorted(beds.items()):
        bs = BedSummary(bb)
        table[(type, "Features")] = bs.nfeats
        table[(type, "Unique bases")] = bs.unique_bases
        table[(type, "Total bases")] = bs.total_bases

    print(tabulate(table), file=sys.stdout)


def gb(args):
    """
    %prog gb gffile fastafile

    Convert GFF3 to Genbank format. Recipe taken from:
    <http://www.biostars.org/p/2492/>
    """
    from Bio.Alphabet import generic_dna

    try:
        from BCBio import GFF
    except ImportError:
        print(
            "You need to install dep first: $ easy_install bcbio-gff", file=sys.stderr
        )

    p = OptionParser(gb.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gff_file, fasta_file = args
    pf = op.splitext(gff_file)[0]
    out_file = pf + ".gb"
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(gff_iter, out_file, "genbank")


def orient(args):
    """
    %prog orient in.gff3 features.fasta > out.gff3

    Change the feature orientations based on translation. This script is often
    needed in fixing the strand information after mapping RNA-seq transcripts.

    You can generate the features.fasta similar to this command:

    $ %prog load --parents=EST_match --children=match_part clc.JCVIv4a.gff
    JCVI.Medtr.v4.fasta -o features.fasta
    """
    from jcvi.formats.fasta import longestorf

    p = OptionParser(orient.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ingff3, fastafile = args
    idsfile = fastafile.rsplit(".", 1)[0] + ".orf.ids"
    if need_update(fastafile, idsfile):
        longestorf([fastafile, "--ids"])

    orientations = DictFile(idsfile)
    gff = Gff(ingff3)
    flipped = 0
    for g in gff:
        id = None
        for tag in ("ID", "Parent"):
            if tag in g.attributes:
                (id,) = g.attributes[tag]
                break
        assert id

        orientation = orientations.get(id, "+")
        if orientation == "-":
            g.strand = {"+": "-", "-": "+"}[g.strand]
            flipped += 1

        print(g)

    logging.debug("A total of {0} features flipped.".format(flipped))


def rename(args):
    """
    %prog rename in.gff3 switch.ids > reindexed.gff3

    Change the IDs within the gff3.
    """
    p = OptionParser(rename.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ingff3, switch = args
    switch = DictFile(switch)

    gff = Gff(ingff3)
    for g in gff:
        (id,) = g.attributes["ID"]
        newname = switch.get(id, id)
        g.attributes["ID"] = [newname]

        if "Parent" in g.attributes:
            parents = g.attributes["Parent"]
            g.attributes["Parent"] = [switch.get(x, x) for x in parents]

        g.update_attributes()
        print(g)


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
        b = next(g.parents(cid, 1))
        print("\t".join((cid, b.id)))


def filter(args):
    """
    %prog filter gffile > filtered.gff

    Filter the gff file based on criteria below:
    (1) feature attribute values: [Identity, Coverage].
    You can get this type of gff by using gmap
    $ gmap -f 2 ....

    (2) Total bp length of child features
    """
    p = OptionParser(filter.__doc__)
    p.add_option(
        "--type", default="mRNA", help="The feature to scan for the attributes"
    )
    g1 = OptionGroup(p, "Filter by identity/coverage attribute values")
    g1.add_option("--id", default=95, type="float", help="Minimum identity")
    g1.add_option("--coverage", default=90, type="float", help="Minimum coverage")
    g1.add_option(
        "--nocase",
        default=False,
        action="store_true",
        help="Case insensitive lookup of attribute names",
    )
    p.add_option_group(g1)
    g2 = OptionGroup(p, "Filter by child feature bp length")
    g2.add_option(
        "--child_ftype", default=None, type="str", help="Child featuretype to consider"
    )
    g2.add_option(
        "--child_bp",
        default=None,
        type="int",
        help="Filter by total bp of children of chosen ftype",
    )
    p.add_option_group(g2)
    p.set_outfile()

    opts, args = p.parse_args(args)
    otype, oid, ocov = opts.type, opts.id, opts.coverage
    cftype, clenbp = opts.child_ftype, opts.child_bp

    id_attr, cov_attr = "Identity", "Coverage"
    if opts.nocase:
        id_attr, cov_attr = id_attr.lower(), cov_attr.lower()

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args

    gffdb = make_index(gffile)
    bad = set()
    ptype = None
    for g in gffdb.features_of_type(otype, order_by=("seqid", "start")):
        if not ptype:
            parent = list(gffdb.parents(g))
            ptype = parent[0].featuretype if len(parent) > 0 else otype
        if cftype and clenbp:
            if gffdb.children_bp(g, child_featuretype=cftype) < clenbp:
                bad.add(g.id)
        elif oid and ocov:
            identity = float(g.attributes[id_attr][0])
            coverage = float(g.attributes[cov_attr][0])
            if identity < oid or coverage < ocov:
                bad.add(g.id)

    logging.debug("{0} bad accns marked.".format(len(bad)))

    fw = must_open(opts.outfile, "w")
    for g in gffdb.features_of_type(ptype, order_by=("seqid", "start")):
        if ptype != otype:
            feats = list(gffdb.children(g, featuretype=otype, order_by="start"))
            ok_feats = [f for f in feats if f.id not in bad]
            if len(ok_feats) > 0:
                g.keep_order = True
                print(g, file=fw)
                for feat in ok_feats:
                    feat.keep_order = True
                    print(feat, file=fw)
                    for child in gffdb.children(feat, order_by="start"):
                        child.keep_order = True
                        print(child, file=fw)
        else:
            if g.id not in bad:
                print(g, file=fw)
                for child in gffdb.children(g, order_by="start"):
                    print(child, file=fw)
    fw.close()


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


def gapsplit(args):
    """
    %prog gapsplit gffile > split.gff

    Read in the gff (normally generated by GMAP) and print it out after splitting
    each feature into one parent and multiple child features based on alignment
    information encoded in CIGAR string.
    """
    p = OptionParser(gapsplit.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args

    gff = Gff(gffile)
    for g in gff:
        if re.match("EST_match", g.type):
            """
            hacky implementation:
            since the standard urlparse.parse_qsl() replaces all "+" symbols with spaces
            we will write a regex to check either for a "-" or a " " (space)
            """
            match = re.search(r"\S+ (\d+) \d+ ([\s{1}\-])", g.attributes["Target"][0])
            if match.group(2) == "-":
                strand = match.group(2)
            else:
                strand = "+"
                g.attributes["Target"][0] = " ".join(
                    str(x) for x in [g.attributes["Target"][0].rstrip(), strand]
                )

            if g.strand == "?":
                g.strand = strand
        else:
            match = re.match(r"\S+ (\d+) \d+", g.attributes["Target"][0])
        target_start = int(match.group(1))

        re_cigar = re.compile(r"(\D+)(\d+)")
        cigar = g.attributes["Gap"][0].split(" ")
        g.attributes["Gap"] = None

        parts = []
        if g.strand == "+":
            for event in cigar:
                match = re_cigar.match(event)
                op, count = match.group(1), int(match.group(2))
                if op in "IHS":
                    target_start += count
                elif op in "DN":
                    g.start += count
                elif op == "P":
                    continue
                else:
                    parts.append(
                        [
                            g.start,
                            g.start + count - 1,
                            target_start,
                            target_start + count - 1,
                        ]
                    )
                    g.start += count
                    target_start += count
        else:
            for event in cigar:
                match = re_cigar.match(event)
                op, count = match.group(1), int(match.group(2))
                if op in "IHS":
                    target_start += count
                elif op in "DN":
                    g.end -= count
                elif op == "P":
                    continue
                else:
                    parts.append(
                        [
                            g.end - count + 1,
                            g.end,
                            target_start,
                            target_start + count - 1,
                        ]
                    )
                    g.end -= count
                    target_start += count

        g.update_attributes()
        print(g)

        parent = g.attributes["Name"][0]
        g.type = "match_part"
        g.attributes.clear()

        for part in parts:
            g.start, g.end = part[0], part[1]
            g.score, g.strand, g.phase = ".", g.strand, "."

            if re.match("EST", g.type):
                target_list = [parent, part[2], part[3], g.strand]
            else:
                target_list = [parent, part[2], part[3]]
            target = " ".join(str(x) for x in target_list)

            g.attributes["Parent"] = [parent]
            g.attributes["Target"] = [target]

            g.update_attributes()
            print(g)


def chain(args):
    """
    %prog chain gffile > chained.gff

    Fill in parent features by chaining child features and return extent of the
    child coordinates.
    """
    valid_merge_op = ("sum", "min", "max", "mean", "collapse")

    p = OptionParser(chain.__doc__)
    p.add_option(
        "--key",
        dest="attrib_key",
        default=None,
        help="Attribute to use as `key` for chaining operation",
    )
    p.add_option(
        "--chain_ftype",
        default="cDNA_match",
        help="GFF feature type to use for chaining operation",
    )
    p.add_option(
        "--parent_ftype",
        default=None,
        help="GFF feature type to use for the chained coordinates",
    )
    p.add_option(
        "--break",
        dest="break_chain",
        action="store_true",
        help="Break long chains which are non-contiguous",
    )
    p.add_option(
        "--transfer_attrib",
        dest="attrib_list",
        help="Attributes to transfer to parent feature; accepts comma"
        + " separated list of attribute names",
    )
    p.add_option(
        "--transfer_score",
        dest="score_merge_op",
        choices=valid_merge_op,
        help="Transfer value stored in score field to parent feature."
        + " Score is reported based on chosen operation",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    attrib_key = opts.attrib_key
    attrib_list = opts.attrib_list
    score_merge_op = opts.score_merge_op
    break_chain = opts.break_chain

    chain_ftype = opts.chain_ftype
    parent_ftype = opts.parent_ftype if opts.parent_ftype else chain_ftype

    gffdict = {}
    fw = must_open(opts.outfile, "w")
    gff = Gff(gffile)
    if break_chain:
        ctr, prev_gid = dict(), None
    for g in gff:
        if g.type != chain_ftype:
            print(g, file=fw)
            continue

        id = g.accn
        gid = id
        if attrib_key:
            assert (
                attrib_key in g.attributes.keys()
            ), "Attribute `{0}` not present in GFF3".format(attrib_key)
            gid = g.get_attr(attrib_key)
        curr_gid = gid
        if break_chain:
            if prev_gid != curr_gid:
                if curr_gid not in ctr:
                    ctr[curr_gid] = 0
                else:
                    ctr[curr_gid] += 1
            gid = "{0}:{1}".format(gid, ctr[curr_gid])
        gkey = (g.seqid, gid)
        if gkey not in gffdict:
            gffdict[gkey] = {
                "seqid": g.seqid,
                "source": g.source,
                "strand": g.strand,
                "type": parent_ftype,
                "coords": [],
                "children": [],
                "score": [],
                "attrs": DefaultOrderedDict(set),
            }
            gffdict[gkey]["attrs"]["ID"].add(gid)

        if attrib_list:
            for a in attrib_list.split(","):
                if a in g.attributes:
                    [gffdict[gkey]["attrs"][a].add(x) for x in g.attributes[a]]
                    del g.attributes[a]

        if break_chain:
            _attrib = "Alias" if attrib_list and ("Name" not in attrib_list) else "Name"
            gffdict[gkey]["attrs"][_attrib].add(curr_gid)

        gffdict[gkey]["coords"].append((g.start, g.end))
        if score_merge_op:
            if is_number(g.score):
                gffdict[gkey]["score"].append(float(g.score))
                g.score = "."

        g.attributes["Parent"] = [gid]
        g.attributes["ID"] = ["{0}-{1}".format(gid, len(gffdict[gkey]["children"]) + 1)]
        g.type = valid_gff_parent_child[g.type]
        g.update_attributes()
        gffdict[gkey]["children"].append(g)
        if break_chain:
            prev_gid = curr_gid

    for gkey, v in sorted(gffdict.items()):
        gseqid, key = gkey
        seqid = v["seqid"]
        source = v["source"]
        type = v["type"]
        strand = v["strand"]
        start, stop = range_minmax(gffdict[gkey]["coords"])

        score = "."
        if score_merge_op:
            v["score"].sort()
            if score_merge_op == "sum":
                score = sum(v["score"])
            elif score_merge_op == "min":
                score = min(v["score"])
            elif score_merge_op == "max":
                score = max(v["score"])
            elif score_merge_op == "mean":
                score = sum(v["score"], 0.0) / len(v["score"])
            elif score_merge_op == "collapse":
                score = ",".join((str(x) for x in v["score"]))

        g = GffLine(
            "\t".join(
                str(x)
                for x in [seqid, source, type, start, stop, score, strand, ".", None]
            )
        )
        g.attributes = v["attrs"]
        g.update_attributes()

        print(g, file=fw)

        for child in gffdict[gkey]["children"]:
            print(child, file=fw)

    fw.close()


def format(args):
    """
    %prog format gffile > formatted.gff

    Read in the gff and print it out, changing seqid, etc.
    """
    from jcvi.formats.obo import GODag_from_SO, validate_term

    valid_multiparent_ops = ["split", "merge"]

    p = OptionParser(format.__doc__)

    g1 = OptionGroup(p, "Parameter(s) used to modify GFF attributes (9th column)")
    g1.add_option("--name", help="Add Name attribute from two-column file")
    g1.add_option("--note", help="Add Note attribute from two-column file")
    g1.add_option(
        "--add_attribute",
        dest="attrib_files",
        help="Add new attribute(s) "
        + "from two-column file(s); attribute name comes from filename; "
        + "accepts comma-separated list of files",
    )
    g1.add_option(
        "--add_dbxref",
        dest="dbxref_files",
        help="Add new Dbxref value(s) (DBTAG:ID) "
        + "from two-column file(s). DBTAG comes from filename, ID comes from 2nd column; "
        + "accepts comma-separated list of files",
    )
    g1.add_option(
        "--nostrict",
        default=False,
        action="store_true",
        help="Disable strict parsing of GFF file and/or mapping file",
    )
    g1.add_option(
        "--remove_attr",
        dest="remove_attrs",
        help="Specify attributes to remove; "
        + "accepts comma-separated list of attribute names",
    )
    g1.add_option(
        "--copy_id_attr_to_name",
        default=False,
        action="store_true",
        help="Copy `ID` attribute value to `Name`, when `Name` is not defined",
    )
    g1.add_option(
        "--invent_name_attr",
        default=False,
        action="store_true",
        help="Invent `Name` attribute for 2nd level child features; "
        + "Formatted like PARENT:FEAT_TYPE:FEAT_INDEX",
    )
    g1.add_option(
        "--no_keep_attr_order",
        default=False,
        action="store_true",
        help="Do not maintain attribute order",
    )
    p.add_option_group(g1)

    g2 = OptionGroup(p, "Parameter(s) used to modify content within columns 1-8")
    g2.add_option(
        "--seqid",
        help="Switch seqid from two-column file. If not"
        + " a file, value will globally replace GFF seqid",
    )
    g2.add_option(
        "--source",
        help="Switch GFF source from two-column file. If not"
        + " a file, value will globally replace GFF source",
    )
    g2.add_option(
        "--type",
        help="Switch GFF feature type from two-column file. If not"
        + " a file, value will globally replace GFF type",
    )
    g2.add_option(
        "--fixphase",
        default=False,
        action="store_true",
        help="Change phase 1<->2, 2<->1",
    )
    p.add_option_group(g2)

    g3 = OptionGroup(
        p, "Other parameter(s) to perform manipulations to the GFF " + "file content"
    )
    g3.add_option(
        "--unique", default=False, action="store_true", help="Make IDs unique"
    )
    g3.add_option(
        "--chaindup",
        default=None,
        dest="duptype",
        help="Chain duplicate features of a particular GFF3 `type`,"
        + " sharing the same ID attribute",
    )
    g3.add_option(
        "--multiparents",
        default=None,
        choices=valid_multiparent_ops,
        help="Split/merge identical features (same `seqid`, `source`, `type`, `coord-range`, `strand`, `phase`) mapping to multiple parents",
    )
    g3.add_option(
        "--remove_feats", help="Comma separated list of features to remove by type"
    )
    g3.add_option(
        "--remove_feats_by_ID",
        help="List of features to remove by ID;"
        + " accepts comma-separated list or list file",
    )
    g3.add_option(
        "--gsac",
        default=False,
        action="store_true",
        help="Fix GSAC GFF3 file attributes",
    )
    g3.add_option(
        "--invent_protein_feat",
        default=False,
        action="store_true",
        help="Invent a protein feature span (chain CDS feats)",
    )
    g3.add_option(
        "--process_ftype",
        default=None,
        type="str",
        help="Specify feature types to process; "
        "accepts comma-separated list of feature types",
    )
    g3.add_option(
        "--gff3", default=False, action="store_true", help="Print output in GFF3 format"
    )
    g3.add_option(
        "--make_gff_store",
        default=False,
        action="store_true",
        help="Store entire GFF file in memory during first iteration",
    )
    p.add_option_group(g3)

    p.set_outfile()
    p.set_SO_opts()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    mapfile = opts.seqid
    names = opts.name
    note = opts.note
    source = opts.source
    ftype = opts.type
    attrib_files = opts.attrib_files.split(",") if opts.attrib_files else None
    dbxref_files = opts.dbxref_files.split(",") if opts.dbxref_files else None
    remove_attrs = opts.remove_attrs.split(",") if opts.remove_attrs else None
    process_ftype = opts.process_ftype.split(",") if opts.process_ftype else None
    gsac = opts.gsac
    assert not (
        opts.unique and opts.duptype
    ), "Cannot use `--unique` and `--chaindup` together"
    assert not (
        opts.type and opts.duptype
    ), "Cannot use `--type` and `--chaindup` together"
    unique = opts.unique
    duptype = opts.duptype
    fixphase = opts.fixphase
    phaseT = {"1": "2", "2": "1"}
    remove_feats = opts.remove_feats.split(",") if opts.remove_feats else None
    remove_feats_by_ID = None
    if opts.remove_feats_by_ID:
        remove_feats_by_ID = (
            LineFile(opts.remove_feats_by_ID, load=True).lines
            if op.isfile(opts.remove_feats_by_ID)
            else opts.remove_feats_by_ID.split(",")
        )
    strict = False if opts.nostrict else True
    make_gff_store = True if gffile in ("-", "stdin") else opts.make_gff_store
    assert not (
        opts.copy_id_attr_to_name and opts.invent_name_attr
    ), "Cannot use `--copy_id_attr_to_name` and `--invent_name_attr` together"
    copy_id_attr_to_name = opts.copy_id_attr_to_name
    invent_name_attr = opts.invent_name_attr
    invent_protein_feat = opts.invent_protein_feat
    compute_signature = False

    outfile = opts.outfile

    mapping = None
    mod_attrs = set()
    if mapfile and op.isfile(mapfile):
        mapping = DictFile(mapfile, delimiter="\t", strict=strict)
        mod_attrs.add("ID")
    if note:
        note = DictFile(note, delimiter="\t", strict=strict)
        mod_attrs.add("Note")
    if source and op.isfile(source):
        source = DictFile(source, delimiter="\t", strict=strict)
    if ftype and op.isfile(ftype):
        ftype = DictFile(ftype, delimiter="\t", strict=strict)
    if names:
        names = DictFile(names, delimiter="\t", strict=strict)
        mod_attrs.add("Name")

    if attrib_files:
        attr_values = {}
        for fn in attrib_files:
            attr_name = op.basename(fn).rsplit(".", 1)[0]
            if attr_name not in reserved_gff_attributes:
                attr_name = attr_name.lower()
            attr_values[attr_name] = DictFile(fn, delimiter="\t", strict=strict)
            mod_attrs.add(attr_name)
    if dbxref_files:
        dbxref_values = {}
        for fn in dbxref_files:
            dbtag = op.basename(fn).rsplit(".", 1)[0]
            dbxref_values[dbtag] = DictFile(fn, delimiter="\t", strict=strict)
        mod_attrs.add("Dbxref")

    if remove_attrs:
        mod_remove_attrs = []
        for remove_attr in remove_attrs:
            if remove_attr in mod_attrs:
                mod_remove_attrs.append(remove_attr)

        if mod_remove_attrs:
            logging.error(
                "Attributes `{0}` cannot be removed and modified".format(
                    ",".join(mod_remove_attrs)
                )
            )
            sys.exit()

    if gsac:  # setting gsac will force IDs to be unique
        unique = True
        notes = {}

    remove = set()
    if (
        unique
        or duptype
        or remove_feats
        or remove_feats_by_ID
        or opts.multiparents == "merge"
        or invent_name_attr
        or make_gff_store
        or invent_protein_feat
    ):
        if unique:
            dupcounts = defaultdict(int)
            seen = defaultdict(int)
            newparentid = {}
        elif duptype:
            dupranges = AutoVivification()
            skip = defaultdict(int)
        if opts.multiparents == "merge":
            merge_feats = AutoVivification()
        if invent_name_attr:
            ft = GffFeatureTracker()
        elif copy_id_attr_to_name:
            pass
        if invent_protein_feat:
            cds_track = {}
        if opts.multiparents == "merge" or invent_name_attr:
            make_gff_store = compute_signature = True
        gff = Gff(
            gffile,
            keep_attr_order=(not opts.no_keep_attr_order),
            make_gff_store=make_gff_store,
            compute_signature=compute_signature,
            strict=strict,
        )
        for g in gff:
            if process_ftype and g.type not in process_ftype:
                continue
            id = g.accn
            if remove_feats and g.type in remove_feats:
                remove.add(id)
            if remove_feats_by_ID and id in remove_feats_by_ID:
                remove.add(id)
            if unique:
                dupcounts[id] += 1
            elif duptype and g.type == duptype:
                dupranges[g.seqid][id][g.idx] = (g.start, g.end)
            if opts.multiparents == "merge" and g.type != "CDS":  # don't merge CDS
                pp = g.get_attr("Parent", first=False)
                if pp and len(pp) > 0:
                    for parent in pp:
                        if parent not in remove:
                            sig = g.sign
                            if sig not in merge_feats:
                                merge_feats[sig]["parents"] = []
                            if parent not in merge_feats[sig]["parents"]:
                                merge_feats[sig]["parents"].append(parent)
            if invent_name_attr:
                parent, iso = atg_name(g.get_attr("Parent"), retval="locus,iso")
                if not parent:
                    parent = g.get_attr("Parent")
                ft.track(parent, g)
            if invent_protein_feat:
                if g.type == "CDS":
                    cds_parent = g.get_attr("Parent")
                    if cds_parent not in cds_track:
                        cds_track[cds_parent] = []
                    cds_track[cds_parent].append((g.start, g.end))

    if opts.verifySO:
        so, _ = GODag_from_SO()
        valid_soterm = {}

    fw = must_open(outfile, "w")
    if not make_gff_store:
        gff = Gff(gffile, keep_attr_order=(not opts.no_keep_attr_order), strict=strict)
    for g in gff:
        if process_ftype and g.type not in process_ftype:
            print(g, file=fw)
            continue

        id = g.accn

        if opts.multiparents == "merge" and g.type != "CDS":  # don't merge CDS
            sig = g.sign
            if len(merge_feats[sig]["parents"]) > 1:
                if "candidate" not in merge_feats[sig]:
                    merge_feats[sig]["candidate"] = id
                    g.set_attr("Parent", merge_feats[sig]["parents"])
                else:
                    continue

        if remove_feats or remove_feats_by_ID:
            if id in remove:
                continue
            else:
                if "Parent" in g.attributes:
                    keep, parent = [], g.get_attr("Parent", first=False)
                    for i, pid in enumerate(parent):
                        if pid not in remove:
                            keep.append(parent[i])
                        else:
                            remove.add(id)
                    if len(keep) == 0:
                        continue
                    parent = g.set_attr("Parent", keep)

        if remove_attrs:
            for remove_attr in remove_attrs:
                if remove_attr in g.attributes:
                    g.set_attr(remove_attr, None)

        if opts.verifySO:
            if g.type not in valid_soterm:
                valid_soterm[g.type] = validate_term(
                    g.type, so=so, method=opts.verifySO
                )
            ntype = valid_soterm[g.type]
            if ntype and g.type != ntype:
                g.type = ntype

        origid = g.seqid
        if fixphase:
            phase = g.phase
            g.phase = phaseT.get(phase, phase)

        if mapfile:
            if isinstance(mapping, dict):
                if origid in mapping:
                    g.seqid = mapping[origid]
                else:
                    logging.error(
                        "{0} not found in `{1}`. ID unchanged.".format(origid, mapfile)
                    )
            else:
                g.seqid = mapfile

        if source:
            if isinstance(source, dict) and g.source in source:
                g.source = source[g.source]
            else:
                g.source = source

        if names:
            if id in names:
                g.set_attr("Name", names[id])

        if note:
            name = g.get_attr("Name")
            tag = None
            if id in note:
                tag = note[id]
            elif name and name in note:
                tag = note[name]

            if tag:
                g.set_attr("Note", tag, update=False)

        if attrib_files:
            for attr_name in attr_values:
                name = g.get_attr("Name")
                if id in attr_values[attr_name]:
                    g.set_attr(attr_name, attr_values[attr_name][id])
                elif name and name in attr_values[attr_name]:
                    g.set_attr(attr_name, attr_values[attr_name][name])

        if dbxref_files:
            for dbtag in dbxref_values:
                if id in dbxref_values[dbtag]:
                    g.set_attr(
                        "Dbxref", dbxref_values[dbtag][id], dbtag=dbtag, append=True
                    )

        if unique:
            if dupcounts[id] > 1:
                seen[id] += 1
                old_id = id
                id = "{0}-{1}".format(old_id, seen[old_id])
                newparentid[old_id] = id
                g.set_attr("ID", id)

            if "Parent" in g.attributes:
                parent = g.attributes["Parent"][0]
                if dupcounts[parent] > 1:
                    g.set_attr("Parent", newparentid[parent])

        if duptype:
            if duptype == g.type and len(dupranges[g.seqid][id]) > 1:
                p = sorted(dupranges[g.seqid][id])
                s, e = dupranges[g.seqid][id][p[0]][
                    0:2
                ]  # get coords of first encountered feature
                if g.start == s and g.end == e and p[0] == g.idx:
                    r = [dupranges[g.seqid][id][x] for x in dupranges[g.seqid][id]]
                    g.start, g.end = range_minmax(r)
                else:
                    skip[(g.seqid, g.idx, id, g.start, g.end)] = 1

        if gsac and g.type == "gene":
            notes[id] = g.attributes["Name"]

        if ftype:
            if isinstance(ftype, dict) and g.type in ftype:
                g.type = ftype[g.type]
            else:
                g.type = ftype

        if invent_name_attr:
            ft.store_symbol(g)
            if re.search(ft.ftype, g.type):
                parent, iso = atg_name(g.get_attr("Parent"), retval="locus,iso")
                if not parent:
                    parent = g.get_attr("Parent")
                if parent in ft.tracker:
                    fidx = ft.feat_index(
                        parent, g.type, g.strand, (g.start, g.end, g.sign)
                    )
                    symbol = ft.get_symbol(parent)
                    attr = "ID" if symbol == parent else "Name"
                    g.set_attr(attr, "{0}:{1}:{2}".format(symbol, g.type, fidx + 1))
                    if opts.multiparents == "merge" and attr == "Name":
                        g.set_attr("ID", "{0}:{1}:{2}".format(parent, g.type, fidx + 1))
        elif copy_id_attr_to_name:
            if "Name" not in g.attributes.keys():
                g.set_attr("Name", g.get_attr("ID"))

        protein_feat = None
        if invent_protein_feat:
            if g.type == "mRNA":
                if id in cds_track:
                    pstart, pstop = range_minmax(cds_track[id])
                    protein_feat = GffLine(
                        "\t".join(
                            str(x)
                            for x in [
                                g.seqid,
                                g.source,
                                "protein",
                                pstart,
                                pstop,
                                ".",
                                g.strand,
                                ".",
                                "ID={0}-Protein;Name={0};Derives_from={0}".format(id),
                            ]
                        )
                    )
            elif g.type == "CDS":
                parent = g.get_attr("Parent")
                if parent in cds_track:
                    _parent = [parent, "{0}-Protein".format(parent)]
                    g.set_attr("Parent", _parent)

        pp = g.get_attr("Parent", first=False)
        if (
            opts.multiparents == "split" and (pp and len(pp) > 1) and g.type != "CDS"
        ):  # separate features with multiple parents
            id = g.get_attr("ID")
            for i, parent in enumerate(pp):
                if id:
                    g.set_attr("ID", "{0}-{1}".format(id, i + 1))
                g.set_attr("Parent", parent, update=True, urlquote=True)
                if gsac:
                    fix_gsac(g, notes)
                print(g, file=fw)
        else:
            if g.gff3 and not opts.gff3:
                opts.gff3 = True
            g.update_attributes(gff3=opts.gff3)
            if gsac:
                fix_gsac(g, notes)
            if duptype == g.type and skip[(g.seqid, g.idx, id, g.start, g.end)] == 1:
                continue
            print(g, file=fw)
            if g.type == "mRNA" and invent_protein_feat:
                print(protein_feat, file=fw)

    fw.close()


def fixboundaries(args):
    """
    %prog fixboundaries gffile --type="gene" --child_ftype="mRNA" > gffile.fixed

    Adjust the boundary coordinates of parents features based on
    range chained child features, extracting their min and max values
    """
    p = OptionParser(fixboundaries.__doc__)
    p.add_option(
        "--type",
        default="gene",
        type="str",
        help="Feature type for which to adjust boundaries",
    )
    p.add_option(
        "--child_ftype",
        default="mRNA",
        type="str",
        help="Child featuretype(s) to use for identifying boundaries",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    gffdb = make_index(gffile)

    fw = must_open(opts.outfile, "w")
    for f in gffdb.all_features(order_by=("seqid", "start")):
        if f.featuretype == opts.type:
            child_coords = []
            for cftype in opts.child_ftype.split(","):
                for c in gffdb.children(f, featuretype=cftype, order_by="start"):
                    child_coords.append((c.start, c.stop))
            f.start, f.stop = range_minmax(child_coords)

        print(f, file=fw)

    fw.close()


def liftover(args):
    """
    %prog liftover gffile > liftover.gff

    Adjust gff coordinates based on tile number. For example,
    "gannotation.asmbl.000095.7" is the 8-th tile on asmbl.000095.
    """
    p = OptionParser(liftover.__doc__)
    p.add_option("--tilesize", default=50000, type="int", help="The size for each tile")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    gff = Gff(gffile)
    for g in gff:
        seqid = g.seqid
        seqid, tilenum = seqid.rsplit(".", 1)
        tilenum = int(tilenum)
        g.seqid = seqid
        offset = tilenum * opts.tilesize
        g.start += offset
        g.end += offset
        print(g)


def get_piles(allgenes):
    """
    Before running uniq, we need to compute all the piles. The piles are a set
    of redundant features we want to get rid of. Input are a list of GffLines
    features. Output are list of list of features distinct "piles".
    """
    from jcvi.utils.range import Range, range_piles

    ranges = [Range(a.seqid, a.start, a.end, 0, i) for i, a in enumerate(allgenes)]

    for pile in range_piles(ranges):
        yield [allgenes[x] for x in pile]


def match_span(f1, f2):
    return (f1.start == f2.start) and (f1.stop == f2.stop)


def match_ftype(f1, f2):
    return f1.featuretype == f2.featuretype


def match_nchildren(f1c, f2c):
    return len(f1c) == len(f2c)


def match_child_ftype(f1c, f2c):
    from collections import Counter

    return len(
        set(Counter(i.featuretype for i in f1c).keys())
        ^ set(Counter(i.featuretype for i in f2c).keys())
    )


def match_Nth_child(f1c, f2c, N=1, slop=False):
    i = N - 1
    f1, f2 = f1c[i], f2c[i]

    if slop:
        if 1 == len(f1c):
            if f1.featuretype.endswith("UTR"):
                if f1.strand == "+":
                    Npos = "F" if f1.featuretype.startswith("five_prime") else "L"
                elif f1.strand == "-":
                    Npos = "L" if f1.featuretype.startswith("five_prime") else "F"
            elif f1.featuretype == "exon":
                return not match_span(f1, f2)
        elif N == 1:
            Npos = "F"
        elif N == len(f1c):
            Npos = "L"

        if Npos == "F":
            return f1.stop == f2.stop
        elif Npos == "L":
            return f1.start == f2.start

    return match_span(f1, f2)


def match_subfeats(f1, f2, dbx1, dbx2, featuretype=None, slop=False):
    """
    Given 2 gffutils features located in 2 separate gffutils databases,
    iterate through all subfeatures of a certain type and check whether
    they are identical or not

    The `slop` parameter allows for variation in the terminal UTR region
    """
    f1c, f2c = (
        list(dbx1.children(f1, featuretype=featuretype, order_by="start")),
        list(dbx2.children(f2, featuretype=featuretype, order_by="start")),
    )

    lf1c, lf2c = len(f1c), len(f2c)
    if match_nchildren(f1c, f2c):
        if lf1c > 0 and lf2c > 0:
            exclN = set()
            if featuretype.endswith("UTR") or featuretype == "exon":
                N = []
                if featuretype.startswith("five_prime"):
                    N = [1] if f1.strand == "+" else [lf1c]
                elif featuretype.startswith("three_prime"):
                    N = [lf1c] if f1.strand == "+" else [1]
                else:  # infer UTR from exon collection
                    N = [1] if 1 == lf1c else [1, lf1c]

                for n in N:
                    if match_Nth_child(f1c, f2c, N=n, slop=slop):
                        exclN.add(n - 1)
                    else:
                        return False

            for i, (cf1, cf2) in enumerate(zip(f1c, f2c)):
                if i in exclN:
                    continue
                if not match_span(cf1, cf2):
                    return False
    else:
        if (lf1c, lf2c) in [(0, 1), (1, 0)] and slop and featuretype.endswith("UTR"):
            return True

        return False

    return True


def import_feats(gffile, type="gene"):
    gff = Gff(gffile)
    allgenes = []
    for g in gff:
        if g.type != type:
            continue
        allgenes.append(g)

    logging.debug("A total of {0} {1} features imported.".format(len(allgenes), type))
    allgenes.sort(key=lambda x: (x.seqid, x.start))
    return allgenes


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
    p.add_option("--type", default="gene", help="Types of features to non-redundify")
    p.add_option("--mode", default="span", choices=supported_modes, help="Pile mode")
    p.add_option("--best", default=1, type="int", help="Use best N features")
    p.add_option(
        "--name",
        default=False,
        action="store_true",
        help="Non-redundify Name attribute",
    )
    p.add_option(
        "--iter",
        default="2",
        choices=("1", "2"),
        help="Number of iterations to grab children",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    mode = opts.mode
    bestn = opts.best

    allgenes = import_feats(gffile, opts.type)
    g = get_piles(allgenes)

    bestids = set()
    for group in g:
        if mode == "span":
            sort_key = lambda x: -x.span
        else:
            sort_key = lambda x: -float(x.score)

        group.sort(key=sort_key)
        seen = set()
        for x in group:
            if len(seen) >= bestn:
                break

            name = x.attributes["Name"][0] if opts.name else x.accn
            if name in seen:
                continue

            seen.add(name)
            bestids.add(x.accn)

    populate_children(opts.outfile, bestids, gffile, iter=opts.iter)


def populate_children(outfile, ids, gffile, iter="2", types=None):
    ids = set(ids)
    fw = must_open(outfile, "w")
    logging.debug("A total of {0} features selected.".format(len(ids)))
    logging.debug("Populate children. Iteration 1..")
    gff = Gff(gffile)
    children = set()
    for g in gff:
        if types and g.type in types:
            ids.add(g.accn)
        if "Parent" not in g.attributes:
            continue
        for parent in g.attributes["Parent"]:
            if parent in ids:
                children.add(g.accn)

    if iter == "2":
        logging.debug("Populate grand children. Iteration 2..")
        gff = Gff(gffile)
        for g in gff:
            if "Parent" not in g.attributes:
                continue
            for parent in g.attributes["Parent"]:
                if parent in children:
                    children.add(g.accn)

    logging.debug("Populate parents..")
    gff = Gff(gffile)
    parents = set()
    for g in gff:
        if g.accn not in ids:
            continue
        if "Parent" not in g.attributes:
            continue
        for parent in g.attributes["Parent"]:
            parents.add(parent)

    combined = ids | children | parents
    logging.debug("Original: {0}".format(len(ids)))
    logging.debug("Children: {0}".format(len(children)))
    logging.debug("Parents: {0}".format(len(parents)))
    logging.debug("Combined: {0}".format(len(combined)))

    logging.debug("Filter gff file..")
    gff = Gff(gffile)
    seen = set()
    for g in gff:
        accn = g.accn
        if accn in seen:
            continue
        if accn in combined:
            seen.add(accn)
            print(g, file=fw)
    fw.close()


def sort(args):
    """
    %prog sort gffile

    Sort gff file using plain old unix sort based on [chromosome, start coordinate].
    or topologically based on hierarchy of features using the gt (genometools) toolkit
    """
    valid_sort_methods = ("unix", "topo")

    p = OptionParser(sort.__doc__)
    p.add_option(
        "--method",
        default="unix",
        choices=valid_sort_methods,
        help="Specify sort method",
    )
    p.add_option(
        "-i",
        dest="inplace",
        default=False,
        action="store_true",
        help="If doing a unix sort, perform sort inplace",
    )
    p.set_tmpdir()
    p.set_outfile()
    p.set_home("gt")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    sortedgff = opts.outfile
    if opts.inplace:
        if opts.method == "topo" or (
            opts.method == "unix" and gffile in ("-", "stdin")
        ):
            logging.error(
                "Cannot perform inplace sort when method is `topo`"
                + " or method is `unix` and input is `stdin` stream"
            )
            sys.exit()

    if opts.method == "unix":
        cmd = "sort"
        cmd += " -k1,1 -k4,4n {0}".format(gffile)
        if opts.tmpdir:
            cmd += " -T {0}".format(opts.tmpdir)
        if opts.inplace:
            cmd += " -o {0}".gffile
            sortedgff = None
        sh(cmd, outfile=sortedgff)
    elif opts.method == "topo":
        GT_HOME = opts.gt_home
        if not op.isdir(GT_HOME):
            logging.error("GT_HOME={0} directory does not exist".format(GT_HOME))
            sys.exit()
        cmd = "{0}".format(op.join(GT_HOME, "bin", "gt"))
        cmd += " gff3 -sort -tidy -retainids -addids no {0}".format(gffile)
        sh(cmd, outfile=sortedgff)


def fromgtf(args):
    """
    %prog fromgtf gtffile

    Convert gtf to gff file. In gtf, the "transcript_id" will convert to "ID=",
    the "transcript_id" in exon/CDS feature will be converted to "Parent=".
    """
    p = OptionParser(fromgtf.__doc__)
    p.add_option(
        "--transcript_id", default="transcript_id", help="Field name for transcript"
    )
    p.add_option("--gene_id", default="gene_id", help="Field name for gene")
    p.add_option(
        "--augustus", default=False, action="store_true", help="Input is AUGUSTUS gtf"
    )
    p.set_home("augustus")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gtffile,) = args
    outfile = opts.outfile
    if opts.augustus:
        ahome = opts.augustus_home
        s = op.join(ahome, "scripts/gtf2gff.pl")
        cmd = "{0} --gff3 < {1} --out={2}".format(s, gtffile, outfile)
        sh(cmd)
        return

    gff = Gff(gtffile)
    fw = must_open(outfile, "w")
    transcript_id = opts.transcript_id
    gene_id = opts.gene_id
    nfeats = 0
    for g in gff:
        if g.type in ("transcript", "mRNA"):
            g.type = "mRNA"
            g.update_tag(transcript_id, "ID")
            g.update_tag("mRNA", "ID")
            g.update_tag(gene_id, "Parent")
            g.update_tag("Gene", "Parent")
        elif g.type in ("exon", "CDS") or "UTR" in g.type:
            g.update_tag("transcript_id", "Parent")
            g.update_tag(g.type, "Parent")
        elif g.type == "gene":
            g.update_tag(gene_id, "ID")
            g.update_tag("Gene", "ID")
        else:
            assert 0, "Don't know how to deal with {0}".format(g.type)

        g.update_attributes()
        print(g, file=fw)
        nfeats += 1

    logging.debug("A total of {0} features written.".format(nfeats))


def frombed(args):
    """
    %prog frombed bed_file [--options] > gff_file

    Convert bed to gff file. In bed, the accn will convert to key='ID'
    Default type will be `match` and default source will be `source`
    """
    p = OptionParser(frombed.__doc__)
    p.add_option("--type", default="match", help="GFF feature type")
    p.add_option("--source", default="default", help="GFF source qualifier")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    bed = Bed(bedfile)

    for b in bed:
        print(b.gffline(type=opts.type, source=opts.source))


def fromsoap(args):
    """
    %prog fromsoap soapfile > gff_file

    """
    p = OptionParser(fromsoap.__doc__)
    p.add_option("--type", default="nucleotide_match", help="GFF feature type")
    p.add_option("--source", default="soap", help="GFF source qualifier")
    p.set_fixchrnames(orgn="maize")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (soapfile,) = args
    pad0 = len(str(sum(1 for line in open(soapfile))))

    fw = must_open(opts.outfile, "w")
    fp = must_open(soapfile)
    for idx, line in enumerate(fp):
        if opts.fix_chr_name:
            from jcvi.utils.cbook import fixChromName

            line = fixChromName(line, orgn=opts.fix_chr_name)

        atoms = line.strip().split("\t")
        attributes = "ID=match{0};Name={1}".format(str(idx).zfill(pad0), atoms[0])
        start, end = int(atoms[8]), int(atoms[5]) + int(atoms[8]) - 1
        seqid = atoms[7]

        print(
            "\t".join(
                str(x)
                for x in (
                    seqid,
                    opts.source,
                    opts.type,
                    start,
                    end,
                    ".",
                    atoms[6],
                    ".",
                    attributes,
                )
            ),
            file=fw,
        )


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

    (gffile,) = args
    gff = Gff(gffile)
    transcript_info = AutoVivification()
    for g in gff:
        if g.type.endswith(("RNA", "transcript")):
            if "ID" in g.attributes and "Parent" in g.attributes:
                transcript_id = g.get_attr("ID")
                gene_id = g.get_attr("Parent")
            elif "mRNA" in g.attributes and "Gene" in g.attributes:
                transcript_id = g.get_attr("mRNA")
                gene_id = g.get_attr("Gene")
            else:
                transcript_id = g.get_attr("ID")
                gene_id = transcript_id
            transcript_info[transcript_id]["gene_id"] = gene_id
            transcript_info[transcript_id]["gene_type"] = g.type
            continue

        if g.type not in valid_gff_to_gtf_type.keys():
            continue

        try:
            transcript_id = g.get_attr("Parent", first=False)
        except IndexError:
            transcript_id = g.get_attr("mRNA", first=False)

        g.type = valid_gff_to_gtf_type[g.type]
        for tid in transcript_id:
            if tid not in transcript_info:
                continue
            gene_type = transcript_info[tid]["gene_type"]
            if not gene_type.endswith("RNA") and not gene_type.endswith("transcript"):
                continue
            gene_id = transcript_info[tid]["gene_id"]
            g.attributes = OrderedDict(
                [("gene_id", [gene_id]), ("transcript_id", [tid])]
            )
            g.update_attributes(gtf=True, urlquote=False)

            print(g)


def merge(args):
    """
    %prog merge gffiles

    Merge several gff files into one. When only one file is given, it is assumed
    to be a file with a list of gff files.
    """
    p = OptionParser(merge.__doc__)
    p.add_option(
        "--seq",
        default=False,
        action="store_true",
        help="Print FASTA sequences at the end",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    nargs = len(args)
    if nargs < 1:
        sys.exit(not p.print_help())

    if nargs == 1:
        (listfile,) = args
        fp = open(listfile)
        gffiles = [x.strip() for x in fp]
    else:
        gffiles = args

    outfile = opts.outfile

    deflines = set()
    fw = must_open(outfile, "w")
    fastarecs = {}
    for gffile in natsorted(gffiles, key=lambda x: op.basename(x)):
        logging.debug(gffile)
        fp = open(gffile)
        for row in fp:
            row = row.rstrip()
            if not row or row[0] == "#":
                if row == FastaTag:
                    break
                if row in deflines:
                    continue
                else:
                    deflines.add(row)

            print(row, file=fw)

        if not opts.seq:
            continue

        f = Fasta(gffile, lazy=True)
        for key, rec in f.iteritems_ordered():
            if key in fastarecs:
                continue
            fastarecs[key] = rec

    if opts.seq:
        print(FastaTag, file=fw)
        SeqIO.write(fastarecs.values(), fw, "fasta")

    fw.close()


def extract(args):
    """
    %prog extract gffile

    --contigs: Extract particular contig(s) from the gff file. If multiple contigs are
    involved, use "," to separate, e.g. "contig_12,contig_150"; or provide a file
    with multiple contig IDs, one per line
    --names: Process particular ID(s) from the gff file. If multiple IDs are
    involved, use "," to separate; or provide a file with multiple IDs, one per line
    """
    p = OptionParser(extract.__doc__)
    p.add_option("--contigs", help="Extract features from certain contigs")
    p.add_option("--names", help="Extract features with certain names")
    p.add_option(
        "--types",
        type="str",
        default=None,
        help="Extract features of certain feature types",
    )
    p.add_option(
        "--children",
        default=0,
        choices=["1", "2"],
        help="Specify number of iterations: `1` grabs children, "
        + "`2` grabs grand-children",
    )
    p.add_option("--tag", default="ID", help="Scan the tags for the names")
    p.add_option(
        "--fasta", default=False, action="store_true", help="Write FASTA if available"
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    contigfile = opts.contigs
    namesfile = opts.names
    typesfile = opts.types
    nametag = opts.tag

    contigID = parse_multi_values(contigfile)
    names = parse_multi_values(namesfile)
    types = parse_multi_values(typesfile)
    outfile = opts.outfile

    if opts.children:
        assert types is not None or names is not None, "Must set --names or --types"
        if names is None:
            names = list()
        populate_children(outfile, names, gffile, iter=opts.children, types=types)
        return

    fp = must_open(gffile)
    fw = must_open(opts.outfile, "w")
    for row in fp:
        atoms = row.split()
        if len(atoms) == 0:
            continue
        tag = atoms[0]
        if row[0] == "#":
            if row.strip() == "###":
                continue
            if not (tag == RegionTag and contigID and atoms[1] not in contigID):
                print(row.rstrip(), file=fw)
            if tag == FastaTag:
                break
            continue

        b = GffLine(row)
        attrib = b.attributes
        if contigID and tag not in contigID:
            continue
        if types and b.type in types:
            _id = b.accn
            if _id not in names:
                names.append(_id)
        if names is not None:
            if nametag not in attrib:
                continue
            if attrib[nametag][0] not in names:
                continue

        print(row.rstrip(), file=fw)

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
        extract([gffile, "--contigs=" + s, "--outfile=" + outfile])


def note(args):
    """
    %prog note gffile > tabfile

    Extract certain attribute field for each feature.
    """
    p = OptionParser(note.__doc__)
    p.add_option(
        "--type",
        default=None,
        help="Only process certain types, multiple types allowed with comma",
    )
    p.add_option(
        "--attribute",
        default="Parent,Note",
        help="Attribute field to extract, multiple fields allowd with comma",
    )
    p.add_option("--AED", type="float", help="Only extract lines with AED score <=")
    p.add_option(
        "--exoncount",
        default=False,
        action="store_true",
        help="Get the exon count for each mRNA feat",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    type = opts.type
    if type:
        type = type.split(",")

    g = make_index(gffile)
    exoncounts = {}
    if opts.exoncount:
        for feat in g.features_of_type("mRNA"):
            nexons = 0
            for c in g.children(feat.id, 1):
                if c.featuretype != "exon":
                    continue
                nexons += 1
            exoncounts[feat.id] = nexons

    attrib = opts.attribute.split(",")

    gff = Gff(gffile)
    seen = set()
    AED = opts.AED
    for g in gff:
        if type and g.type not in type:
            continue
        if AED is not None and float(g.attributes["_AED"][0]) > AED:
            continue
        keyval = [g.accn] + [
            ",".join(g.attributes[x]) for x in attrib if x in g.attributes
        ]
        if exoncounts:
            nexons = exoncounts.get(g.accn, 0)
            keyval.append(str(nexons))
        keyval = tuple(keyval)
        if keyval not in seen:
            print("\t".join(keyval))
            seen.add(keyval)


def splicecov(args):
    """
    %prog splicecov annotation.gff3 junctions.bed

    Given an annotation GFF file (containing introns) and a
    TopHat junctions.bed file (preprocessed using formats.bed.juncs(),
    each intron gets tagged with the JUNC identifier and read coverage.

    Output is a summary table listing for each gene locus, the isoform number,
    number of splice junctions and {average, median, min & max} read coverage
    across the junctions.
    """
    from tempfile import mkstemp
    from pybedtools import BedTool
    from jcvi.utils.cbook import SummaryStats

    p = OptionParser(splicecov.__doc__)
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    (
        gfffile,
        juncsbed,
    ) = args
    tagged = "{0}.{1}.gff3".format(gfffile.rsplit(".", 1)[0], "tag_introns")

    gff3, junc = BedTool(gfffile), BedTool(juncsbed)
    ab = gff3.intersect(junc, wao=True, f=1.0, s=True)
    abfh = must_open(ab.fn)

    seen = set()
    scov = AutoVivification()

    fh, tmpgff = mkstemp(suffix=".gff3")
    fw = must_open(tmpgff, "w")
    for line in abfh:
        args = line.strip().split("\t")
        g = GffLine("\t".join(str(x) for x in args[:9]))
        if g.type == "intron" and args[10] != -1:
            ispan, jspan = g.span, int(args[11]) - int(args[10])
            if ispan == jspan:
                g.set_attr("ID", args[12], update=True)
                g.score = int(args[13])

                pparts = g.get_attr("Parent").split(".")
                locus, iso = pparts[0], ".".join(pparts[1:])
                seen.add(iso)
                if not scov[locus][iso]:
                    scov[locus][iso] = []
                scov[locus][iso].append(g.score)
            else:
                continue
        print(g, file=fw)
    fw.close()

    format([tmpgff, "--unique", "-o", tagged])
    os.unlink(tmpgff)

    isos = sorted(list(seen))
    fw = must_open(opts.outfile, "w")
    h1, h2, stats = ["#"], ["#locus"], ["N", "mean", "median", "min", "max"]
    for iso in isos:
        h1.extend([str(iso)] + [""] * (len(stats) - 1))
        h2.extend(stats)
    print("\t".join(str(x) for x in h1), file=fw)
    print("\t".join(str(x) for x in h2), file=fw)
    for locus in scov.keys():
        out = [locus]
        for iso in isos:
            if iso in scov[locus].keys():
                juncs = scov[locus][iso]
                jstats = SummaryStats(juncs, dtype="int")
                out.extend(
                    [jstats.size, jstats.mean, jstats.median, jstats.min, jstats.max]
                )
            else:
                out.extend(["-"] * len(stats))
        print("\t".join(str(x) for x in out), file=fw)
    fw.close()


def bed(args):
    """
    %prog bed gff_file [--options]

    Parses the start, stop locations of the selected features out of GFF and
    generate a bed file
    """
    p = OptionParser(bed.__doc__)
    p.add_option(
        "--type",
        dest="type",
        default="gene",
        help="Feature type to extract, use comma for multiple",
    )
    p.add_option("--key", default="ID", help="Key in the attributes to extract")
    p.add_option("--accn", help="Use fixed accn in the 4th column")
    p.add_option("--source", help="Source to extract from, use comma for multiple")
    p.add_option(
        "--span",
        default=False,
        action="store_true",
        help="Use feature span in the score column",
    )
    p.add_option(
        "--score_attrib",
        dest="score_attrib",
        default=False,
        help="Attribute whose value is to be used as score in `bedline`",
    )
    p.add_option(
        "--append_source",
        default=False,
        action="store_true",
        help="Append GFF source name to extracted key value",
    )
    p.add_option(
        "--append_ftype",
        default=False,
        action="store_true",
        help="Append GFF feature type to extracted key value",
    )
    p.add_option(
        "--nosort",
        default=False,
        action="store_true",
        help="Do not sort the output bed file",
    )
    p.add_option(
        "--primary_only",
        default=False,
        action="store_true",
        help="Only retains a single transcript per gene",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    key = opts.key or None
    type = opts.type or set()
    accn = opts.accn
    source = opts.source or set()
    span = opts.span
    primary_only = opts.primary_only

    if opts.type:
        type = set(x.strip() for x in opts.type.split(","))
    if opts.source:
        source = set(x.strip() for x in opts.source.split(","))

    gff = Gff(
        gffile,
        key=key,
        append_source=opts.append_source,
        append_ftype=opts.append_ftype,
        score_attrib=opts.score_attrib,
    )
    b = Bed()
    seen_parents = set()  # used with --primary_only
    skipped_non_primary = 0

    for g in gff:
        if type and g.type not in type:
            continue
        if source and g.source not in source:
            continue
        if primary_only:
            if g.parent in seen_parents:
                skipped_non_primary += 1
                continue
            elif g.parent:
                seen_parents.add(g.parent)

        bl = g.bedline
        if accn:
            bl.accn = accn
        if span:
            bl.score = bl.span
        b.append(bl)

    sorted = not opts.nosort
    b.print_to_file(opts.outfile, sorted=sorted)
    logging.debug(
        "Extracted {0} features (type={1} id={2})".format(len(b), ",".join(type), key)
    )
    if primary_only:
        logging.debug("Skipped non-primary: %d", skipped_non_primary)


def make_index(gff_file):
    """
    Make a sqlite database for fast retrieval of features.
    """
    import gffutils

    db_file = gff_file + ".db"

    if need_update(gff_file, db_file):
        if op.exists(db_file):
            os.remove(db_file)
        logging.debug("Indexing `{0}`".format(gff_file))
        gffutils.create_db(gff_file, db_file, merge_strategy="create_unique")
    else:
        logging.debug("Load index `{0}`".format(gff_file))

    return gffutils.FeatureDB(db_file)


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
    p.add_option(
        "--parents",
        default="gene",
        help="list of features to extract, use comma to separate (e.g. 'gene,mRNA')",
    )

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gff_file,) = args
    g = make_index(gff_file)
    parents = set(opts.parents.split(","))

    for feat in get_parents(gff_file, parents):

        cc = [c.id for c in g.children(feat.id, 1)]
        if len(cc) <= 1:
            continue

        print("\t".join(str(x) for x in (feat.id, feat.start, feat.stop, "|".join(cc))))


def load(args):
    """
    %prog load gff_file fasta_file [--options]

    Parses the selected features out of GFF, with subfeatures concatenated.
    For example, to get the CDS sequences, do this:
    $ %prog load athaliana.gff athaliana.fa --parents mRNA --children CDS

    To get 500bp upstream of a genes Transcription Start Site (TSS), do this:
    $ %prog load athaliana.gff athaliana.fa --feature=upstream:TSS:500

    Switch TSS with TrSS for Translation Start Site.
    """
    from datetime import datetime as dt
    from jcvi.formats.fasta import Seq, SeqRecord

    # can request output fasta sequence id to be picked from following attributes
    valid_id_attributes = ["ID", "Name", "Parent", "Alias", "Target", "orig_protein_id"]

    p = OptionParser(load.__doc__)
    p.add_option(
        "--parents",
        dest="parents",
        default="mRNA",
        help="list of features to extract, use comma to separate (e.g."
        + "'gene,mRNA')",
    )
    p.add_option(
        "--children",
        dest="children",
        default="CDS",
        help="list of features to extract, use comma to separate (e.g."
        + "'five_prime_UTR,CDS,three_prime_UTR')",
    )
    p.add_option(
        "--feature",
        dest="feature",
        help="feature type to extract. e.g. `--feature=CDS` or "
        + "`--feature=upstream:TSS:500`",
    )
    p.add_option(
        "--id_attribute",
        choices=valid_id_attributes,
        help="The attribute field to extract and use as FASTA sequence ID",
    )
    p.add_option(
        "--desc_attribute",
        default="Note",
        help="The attribute field to extract and use as FASTA sequence description",
    )
    p.add_option(
        "--full_header",
        default=None,
        choices=["default", "tair"],
        help="Specify if full FASTA header (with seqid, coordinates and datestamp) should be generated",
    )

    g1 = OptionGroup(p, "Optional parameters (if generating full header)")
    g1.add_option(
        "--sep",
        dest="sep",
        default=" ",
        help="Specify separator used to delimiter header elements",
    )
    g1.add_option(
        "--datestamp",
        dest="datestamp",
        help="Specify a datestamp in the format YYYYMMDD or automatically pick `today`",
    )
    g1.add_option(
        "--conf_class",
        dest="conf_class",
        default=False,
        action="store_true",
        help="Specify if `conf_class` attribute should be parsed and placed in the header",
    )
    p.add_option_group(g1)

    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    gff_file, fasta_file = args

    if opts.feature:
        (
            opts.feature,
            opts.parent,
            opts.children,
            upstream_site,
            upstream_len,
            flag,
            error_msg,
        ) = parse_feature_param(opts.feature)
        if flag:
            sys.exit(error_msg)

    parents = set(opts.parents.split(","))
    children_list = set(opts.children.split(","))

    """
    In a situation where we want to extract sequence for only the top-level
    parent feature, specify feature type of parent == child
    """
    skipChildren = (
        True if len(parents.symmetric_difference(children_list)) == 0 else False
    )

    id_attr = opts.id_attribute
    desc_attr = opts.desc_attribute
    sep = opts.sep

    import gffutils

    g = make_index(gff_file)
    f = Fasta(fasta_file, index=False)
    seqlen = {}
    for seqid, size in f.itersizes():
        seqlen[seqid] = size

    fw = must_open(opts.outfile, "w")

    for feat in get_parents(gff_file, parents):
        desc = ""
        if desc_attr:
            fparent = (
                feat.attributes["Parent"][0] if "Parent" in feat.attributes else None
            )
            if fparent:
                try:
                    g_fparent = g[fparent]
                except gffutils.exceptions.FeatureNotFoundError:
                    logging.error("{} not found in index .. skipped".format(fparent))
                    continue
                if desc_attr in g_fparent.attributes:
                    desc = ",".join(g_fparent.attributes[desc_attr])
            if not desc and desc_attr in feat.attributes:
                desc = ",".join(feat.attributes[desc_attr])

        if opts.full_header:
            desc_parts = []
            desc_parts.append(desc)

            if opts.conf_class and "conf_class" in feat.attributes:
                desc_parts.append(feat.attributes["conf_class"][0])

            if opts.full_header == "tair":
                orient = "REVERSE" if feat.strand == "-" else "FORWARD"
                feat_coords = "{0}:{1}-{2} {3} LENGTH=[LEN]".format(
                    feat.seqid, feat.start, feat.end, orient
                )
            else:
                (s, e) = (
                    (feat.start, feat.end)
                    if (feat.strand == "+")
                    else (feat.end, feat.start)
                )
                feat_coords = "{0}:{1}-{2}".format(feat.seqid, s, e)
            desc_parts.append(feat_coords)

            datestamp = (
                opts.datestamp
                if opts.datestamp
                else "{0}{1}{2}".format(dt.now().year, dt.now().month, dt.now().day)
            )
            desc_parts.append(datestamp)

            desc = sep.join(str(x) for x in desc_parts)
            desc = "".join(str(x) for x in (sep, desc)).strip()

        if opts.feature == "upstream":
            upstream_start, upstream_stop = get_upstream_coords(
                upstream_site, upstream_len, seqlen[feat.seqid], feat, children_list, g
            )

            if not upstream_start or not upstream_stop:
                continue

            feat_seq = f.sequence(
                dict(
                    chr=feat.seqid,
                    start=upstream_start,
                    stop=upstream_stop,
                    strand=feat.strand,
                )
            )

            (s, e) = (
                (upstream_start, upstream_stop)
                if feat.strand == "+"
                else (upstream_stop, upstream_start)
            )
            upstream_seq_loc = str(feat.seqid) + ":" + str(s) + "-" + str(e)
            desc = sep.join(
                str(x)
                for x in (desc, upstream_seq_loc, "FLANKLEN=" + str(upstream_len))
            )
        else:
            children = []
            if not skipChildren:
                for c in g.children(feat.id, 1):
                    if c.featuretype not in children_list:
                        continue
                    child = f.sequence(
                        dict(chr=c.chrom, start=c.start, stop=c.stop, strand=c.strand)
                    )
                    children.append((child, c))

                if not children:
                    print(
                        "[warning] %s has no children with type %s"
                        % (feat.id, ",".join(children_list)),
                        file=sys.stderr,
                    )
                    continue
            else:
                child = f.sequence(
                    dict(
                        chr=feat.seqid,
                        start=feat.start,
                        stop=feat.end,
                        strand=feat.strand,
                    )
                )
                children.append((child, feat))

            # sort children in incremental position
            children.sort(key=lambda x: x[1].start)
            # reverse children if negative strand
            if feat.strand == "-":
                children.reverse()
            feat_seq = "".join(x[0] for x in children)

        desc = desc.replace('"', "")

        id = (
            ",".join(feat.attributes[id_attr])
            if id_attr and feat.attributes[id_attr]
            else feat.id
        )

        if opts.full_header == "tair":
            desc = desc.replace("[LEN]", str(len(feat_seq)))

        rec = SeqRecord(Seq(feat_seq), id=id, description=desc)
        SeqIO.write([rec], fw, "fasta")
        fw.flush()


def parse_feature_param(feature):
    """
    Take the --feature param (coming from gff.load() and parse it.
    Returns feature, parents and children terms.

    Also returns length of upstream sequence (and start site) requested

    If erroneous, returns a flag and error message to be displayed on exit
    """
    # can request upstream sequence only from the following valid sites
    valid_upstream_sites = ["TSS", "TrSS"]

    upstream_site, upstream_len = None, None
    flag, error_msg = None, None
    parents, children = None, None
    if re.match(r"upstream", feature):
        parents, children = "mRNA", "CDS"
        feature, upstream_site, upstream_len = re.search(
            r"([A-z]+):([A-z]+):(\S+)", feature
        ).groups()

        if not is_number(upstream_len):
            flag, error_msg = (
                1,
                "Error: upstream len `" + upstream_len + "` should be an integer",
            )

        upstream_len = int(upstream_len)
        if upstream_len < 0:
            flag, error_msg = (
                1,
                "Error: upstream len `" + str(upstream_len) + "` should be > 0",
            )

        if upstream_site not in valid_upstream_sites:
            flag, error_msg = (
                1,
                f"Error: upstream site `{upstream_site}` not valid. Please choose from {valid_upstream_sites}",
            )
    elif feature == "CDS":
        parents, children = "mRNA", "CDS"
    else:
        flag, error_msg = 1, "Error: unrecognized option --feature=" + feature

    return feature, parents, children, upstream_site, upstream_len, flag, error_msg


def get_upstream_coords(uSite, uLen, seqlen, feat, children_list, gffdb):
    """
    Subroutine takes upstream site, length, reference sequence length,
    parent mRNA feature (GffLine object), list of child feature types
    and a GFFutils.GFFDB object as the input

    If upstream of TSS is requested, use the parent feature coords
    to extract the upstream sequence

    If upstream of TrSS is requested,  iterates through all the
    children (CDS features stored in the sqlite GFFDB) and use child
    feature coords to extract the upstream sequence

    If success, returns the upstream start and stop coordinates
    else, returns None
    """
    if uSite == "TSS":
        (upstream_start, upstream_stop) = (
            (feat.start - uLen, feat.start - 1)
            if feat.strand == "+"
            else (feat.end + 1, feat.end + uLen)
        )
    elif uSite == "TrSS":
        children = []
        for c in gffdb.children(feat.id, 1):

            if c.featuretype not in children_list:
                continue
            children.append((c.start, c.stop))

        if not children:
            print(
                "[warning] %s has no children with type %s"
                % (feat.id, ",".join(children_list)),
                file=sys.stderr,
            )
            return None, None

        cds_start, cds_stop = range_minmax(children)
        (upstream_start, upstream_stop) = (
            (cds_start - uLen, cds_start - 1)
            if feat.strand == "+"
            else (cds_stop + 1, cds_stop + uLen)
        )

    if feat.strand == "+" and upstream_start < 1:
        upstream_start = 1
    elif feat.strand == "-" and upstream_stop > seqlen:
        upstream_stop = seqlen

    actual_uLen = upstream_stop - upstream_start + 1
    if actual_uLen < uLen:
        print(
            "[warning] sequence upstream of {0} ({1} bp) is less than upstream length {2}".format(
                feat.id, actual_uLen, uLen
            ),
            file=sys.stderr,
        )
        return None, None

    return upstream_start, upstream_stop


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
    p.add_option("--parent", default="mRNA", help="Top feature type")
    p.add_option("--block", default="exon", help="Feature type for regular blocks")
    p.add_option("--thick", default="CDS", help="Feature type for thick blocks")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
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
        # When there is no thick part, thickStart and thickEnd are usually set
        # to the chromStart position
        # <https://genome.ucsc.edu/FAQ/FAQformat.html>
        thickStart = chromStart
        thickEnd = chromStart
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

        print(
            "\t".join(
                str(x)
                for x in (
                    chrom,
                    chromStart,
                    chromEnd,
                    name,
                    score,
                    strand,
                    thickStart,
                    thickEnd,
                    itemRgb,
                    blockCount,
                    blockSizes,
                    blockStarts,
                )
            ),
            file=fw,
        )


if __name__ == "__main__":
    main()
