#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Genbank AGP file format, see spec here
http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml
"""

import sys
import logging

from copy import deepcopy
from optparse import OptionParser
from collections import defaultdict
from itertools import groupby

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from jcvi.formats.base import LineFile
from jcvi.formats.fasta import Fasta
from jcvi.assembly.base import A50 
from jcvi.utils.iter import pairwise
from jcvi.apps.base import ActionDispatcher, set_debug


Valid_component_type = "ADFGNOPUW"
Valid_gap_type = ("fragment", "clone", "contig", "centromere", "short_arm",
        "heterochromatin", "telomere", "repeat")
Valid_orientation = ("+", "-", "0", "na")


class AGPLine (object):

    def __init__(self, row, validate=True):

        atoms = row.split('\t')
        atoms[-1] = atoms[-1].strip()
        self.object = atoms[0]
        self.object_beg = int(atoms[1])
        self.object_end = int(atoms[2])
        self.object_span = self.object_end - self.object_beg + 1
        self.part_number = atoms[3]
        self.component_type = atoms[4]
        self.is_gap = (self.component_type=='N')
        
        if not self.is_gap:
            self.component_id = atoms[5]
            self.component_beg = int(atoms[6])
            self.component_end = int(atoms[7])
            self.component_span = self.component_end - self.component_beg + 1
            self.orientation = atoms[8].strip()
        else:
            self.gap_length = int(atoms[5])
            self.gap_type = atoms[6]
            self.linkage = atoms[7]
            self.empty = ""
            self.orientation = "na"

        if validate:
            try:
                self.validate()
            except AssertionError, b:
                logging.error("%s\nerror when validating this line:\n%s" % \
                        (b, row))

    def __str__(self):
        
        fields = [self.object, self.object_beg, self.object_end,
                self.part_number, self.component_type]

        if not self.is_gap:
            fields += [self.component_id, self.component_beg,
                    self.component_end, self.orientation]
        else:
            fields += [self.gap_length, self.gap_type, 
                    self.linkage, self.empty]

        return "\t".join(str(x) for x in fields)

    __repr__ = __str__

    @property
    def bedline(self):
        # bed formatted line
        gid = self.component_id if not self.is_gap else self.gap_type 
        return "\t".join((self.object, str(self.object_beg-1),
                str(self.object_end), gid, self.component_type, self.orientation))

    def validate(self):
        assert self.component_type in Valid_component_type, \
                "component_type has to be one of %s" % Valid_component_type
        assert self.object_beg <= self.object_end, \
                "object_beg needs to be <= object_end"

        if not self.is_gap:
            assert self.component_beg <= self.component_end, \
                    "component_begin needs to be <= component_end"
            assert self.object_span==self.component_span, \
                    "object_span (%d) needs to be same as component_span (%d)" %\
                    (self.object_span, self.component_span)
        else:
            assert self.gap_length >= 1, \
                    "gap_length needs to be >= 1"
            assert self.object_span == self.gap_length, \
                    "object span (%d) needs to be same as gap_length (%d)" % \
                    (self.object_span, self.gap_length)


class AGP (LineFile):

    def __init__(self, filename, validate=True):
        super(AGP, self).__init__(filename)

        fp = open(filename)
        for row in fp:
            if row[0]=='#': continue
            self.append(AGPLine(row, validate=validate))

        self.validate = validate
        if validate:
            self.sort(key=lambda x: (x.object, x.object_beg))
            self.validate_all()

    @property
    def simple_agp(self):
        """
        returns a dict with component_id => agpline
        """
        return dict((x.component_id, x) for x in self if not x.is_gap)


    @classmethod
    def print_header(cls, fw=sys.stdout, organism="Medicago truncatula", taxid=3880,
            source="J. Craig Venter Institute", comment=None):
        # these comments are entirely optional, modeled after maize AGP
        print >> fw, "# ORGANISM: {0}".format(organism)
        print >> fw, "# TAX_ID: {0}".format(taxid)
        print >> fw, "# GENOME CENTER: {0}".format(source)
        if comment:
            print >> fw, "# COMMENT: {0}".format(comment)
        header = "object object_beg object_end part_number component_type " +\
                 "component_id/gap_length component_beg/gap_type " +\
                 "component_end/linkage orientation"
        print >> fw, "# FIELDS: {0}".format(", ".join(header.split()))


    def report_stats(self, object, bacs, scaffold_sizes):
        
        from jcvi.utils.cbook import human_size

        nbacs = len(bacs)
        nscaffolds = len(scaffold_sizes)
        a50, l50, n50, scaffold_sizes = A50(scaffold_sizes)

        print "\t".join(str(x) for x in (object, nbacs, nscaffolds, l50,
            human_size(n50, precision=2, target="Mb")))


    def summary_one(self, object, lines):
        bacs = set()
        scaffold_sizes = []
        _scaffold_key = lambda x: x.is_gap and \
                x.gap_type in ("clone", "contig") and x.linkage=="no"

        for is_gap, scaffold in groupby(lines, key=_scaffold_key):
            if is_gap: continue

            scaffold = list(scaffold)
            scaffold_size = 0
            for b in scaffold:
                if b.is_gap:
                    scaffold_size += b.gap_length 
                else:
                    bacs.add(b.component_id)
                    scaffold_size += b.component_span

            scaffold_sizes.append(scaffold_size)

        self.report_stats(object, bacs, scaffold_sizes)

        return bacs, scaffold_sizes


    def summary_all(self):
        
        all_bacs = set()
        all_scaffold_sizes = []
        for ob, lines_with_same_ob in groupby(self, key=lambda x: x.object):
            lines = list(lines_with_same_ob)
            bacs, scaffold_sizes = self.summary_one(ob, lines)
            all_bacs |= bacs
            all_scaffold_sizes.extend(scaffold_sizes)

        self.report_stats("Total", all_bacs, all_scaffold_sizes)


    def validate_one(self, object, lines):
        object_beg = lines[0].object_beg
        assert object_beg==1, \
                "object %s needs to start at 1 (instead of %d)" % \
                (object, object_beg)

        for a, b in pairwise(lines):
            assert b.object_beg-a.object_end==1, \
                    "lines not continuous coords between:\n%s\n%s" % \
                    (a, b)

    def validate_all(self):
        for ob, lines_with_same_ob in groupby(self, key=lambda x: x.object):
            lines = list(lines_with_same_ob)
            self.validate_one(ob, lines)


    def build_one(self, object, lines, fasta, fw):
        """
        Construct molecule using component fasta sequence 
        """
        components = []

        total_bp = 0
        for line in lines:

            if line.is_gap:
                seq = 'N' * line.gap_length
            else:
                seq = fasta.sequence(dict(chr=line.component_id, 
                        start=line.component_beg, 
                        stop=line.component_end,
                        strand=line.orientation))

            components.append(seq)
            total_bp += len(seq) 

            if self.validate:
                assert total_bp == line.object_end, \
                        "cumulative base pairs (%d) does not match (%d)" % \
                        (total_bp, line.object_end)


        rec = SeqRecord(Seq(''.join(components)), id=object, description="")
        SeqIO.write([rec], fw, "fasta")
        logging.debug("Write object %s to `%s`" % (object, fw.name))


    def build_all(self, componentfasta, targetfasta):

        f = Fasta(componentfasta, index=False)
        fw = open(targetfasta, "w")

        for ob, lines_with_same_ob in groupby(self, key=lambda x: x.object):

            lines = list(lines_with_same_ob)
            self.build_one(ob, lines, f, fw)


def main():

    actions = (
        ('summary', 'print out a table of scaffold statistics'),
        ('phase', 'given genbank file, get the phase for the HTG BAC record'),
        ('bed', 'print out the tiling paths in bed format'),
        ('gaps', 'print out the distribution of gap sizes'),
        ('accessions', 'print out a list of accessions'),
        ('chr0', 'build AGP file for unassembled sequences'),
        ('reindex', 'assume the component line order is correct, modify coordinates'),
        ('build', 'given agp file and component fasta file, build the ' +\
                 'pseudomolecule fasta'),
        ('validate', 'given agp file, component fasta and pseudomolecule fasta, ' +\
                     'validate if the build is correct')
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def reindex(args):
    """
    %prog agpfile newagpfile

    assume the component line order is correct, modify coordinates, this is
    necessary mostly due to manual edits (insert/delete component) that disrupts
    the target coordinates
    """
    p = OptionParser(reindex.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    argfile, newagpfile = args
    agp = AGP(argfile, validate=False)
    fw = open(newagpfile, "w")
    for chr, chr_agp in groupby(agp, lambda x: x.object):
        chr_agp = list(chr_agp)
        object_beg = 1
        for i, b in enumerate(chr_agp):
            b.object_beg = object_beg
            b.part_number = i + 1

            if not b.is_gap:
                b.object_end = object_beg + b.component_span - 1
            else:
                b.object_end = object_beg + b.gap_length - 1
            
            object_beg = b.object_end + 1

            print >> fw, str(b)

    # Last step: validate the new agpfile
    fw.close()
    agp = AGP(newagpfile, validate=True)


def summary(args):
    """
    %prog agpfile

    print a table of scaffold statistics, number of BACs, no of scaffolds,
    scaffold N50, scaffold L50, actual sequence, PSMOL NNNs, PSMOL-length, % of
    PSMOL seq'd
    """
    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    agpfile, = args
    header = "Chromosome|Number of BACs|No of Scaffolds|Scaff N50|Scaff L50"
    header = header.replace('|', '\t')
    print header

    agp = AGP(agpfile)
    agp.summary_all()


def get_phase(rec):
    keywords = rec.annotations["keywords"]
    description = rec.description.upper()

    if "HTGS_PHASE1" in keywords:
        phase = 1
    elif "HTGS_PHASE2" in keywords:
        phase = 2
    elif len(keywords)==1 and "HTG" in keywords:
        phase = 3
    elif "PLN" in keywords: # EMBL BACs
        if "DRAFT" in description:
            if "UNORDERED" in description:
                phase = 1
            else:
                phase = 2
        else:
            assert "COMPLETE" in description, description
            phase = 3
    else:
        logging.error("{0}: {1}".format(rec.name, description))
        phase = 4

    return phase, keywords


def phase(args):
    """
    %prog phase genbankfile

    Input has to be gb file. Search the `KEYWORDS` section to look for PHASE.
    """
    p = OptionParser(phase.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    gbfile, = args
    for rec in SeqIO.parse(gbfile, "gb"):
        bac_phase, keywords = get_phase(rec)
        print "{0}\t{1}\t{2}".format(rec.id, bac_phase, "; ".join(keywords))


# phase 0 - (P)refinish; phase 1,2 - (D)raft; phase 3 - (F)inished; 4 - (O)thers
Phases = "PDDFO"

def iter_phase(phasefile):
    fp = open(phasefile)
    for row in fp:
        id, phase, keywords = row.split("\t")
        yield id, Phases[int(phase)]


def chr0(args):
    """
    %prog fastafile [phasefile]

    build AGP file for unassembled sequences, and add gaps between. Phase list
    contains two columns - BAC and phase (0, 1, 2, 3).
    """
    p = OptionParser(chr0.__doc__)
    p.add_option("--gapsize", dest="gapsize", default=0, type="int",
            help="create a new molecule chr0 with x N's inserted between " +\
                 "[default: do not create new molecule]")
    opts, args = p.parse_args(args)

    nargs = len(args)
    if nargs not in (1, 2):
        sys.exit(p.print_help())

    if nargs == 2:
        fastafile, phasefile = args
        phases = dict(iter_phase(phasefile))
    else:
        fastafile, = args
        f = Fasta(fastafile)
        phases = dict((x, 'W') for x in f.iterkeys())

    agpfile = fastafile.rsplit(".", 1)[0] + ".agp"
    f = Fasta(fastafile)
    fw = open(agpfile, "w")

    AGP.print_header(fw, 
        comment="{} components with undetermined chromosomal locations".format(len(f)))

    gap_length = opts.gapsize
    object_beg = 1

    if gap_length:
        object = "chr0"
        gap_type = "clone"
        linkage = "no"

        part_number = 0
        for component_id, size in f.itersizes_ordered():
            if part_number > 0: # print gap except for the first one
                object_end = object_beg + gap_length - 1
                part_number += 1
                print >> fw, "\t".join(str(x) for x in \
                        (object, object_beg, object_end, part_number,
                         'N', gap_length, gap_type, linkage, ""))

                object_beg += gap_length

            object_end = object_beg + size -1
            part_number += 1
            print >> fw, "\t".join(str(x) for x in \
                    (object, object_beg, object_end, part_number,
                     phases[component_id], component_id, 1, size, '0'))

            object_beg += size
    else:

        part_number = 1
        scaffold_number = 0
        for component_id, size in f.itersizes_ordered():
            #object_id = component_id.rsplit(".")[0]
            scaffold_number += 1
            object_id = "scaffold{0:03d}".format(scaffold_number)
            object_end = size
            print >> fw, "\t".join(str(x) for x in \
                    (object_id, object_beg, object_end, part_number,
                    phases[component_id], component_id, 1, size, '0'))


def accessions(args):
    """
    %prog accessions agpfile

    print out a list of accessions, one per line
    """
    p = OptionParser(accessions.__doc__)
    p.add_option("--noversion", dest="noversion", default=False, action="store_true",
            help="Remove trailing accession versions")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    agpfile, = args
    agp = AGP(agpfile)
    seen = set()
    for a in agp:
        if a.is_gap: continue
        component_id = a.component_id

        if opts.noversion:
            component_id = component_id.rsplit(".", 1)[0]

        if component_id not in seen:
            print component_id
        seen.add(component_id)


def bed(args):
    """
    %prog bed agpfile

    print out the tiling paths in bed format
    """
    p = OptionParser(bed.__doc__)
    p.add_option("--nogaps", dest="nogaps", default=False, action="store_true",
            help="do not print bed lines for gaps")
    opts, args = p.parse_args(args)
    
    if len(args) != 1:
        sys.exit(p.print_help())

    agpfile, = args
    agp = AGP(agpfile)
    for a in agp:
        if opts.nogaps and a.is_gap: continue
        print a.bedline


def gaps(args):
    """
    %prog gaps agpfile

    print out the distribution of gapsizes
    """
    p = OptionParser(gaps.__doc__)
    p.add_option("--merge", dest="merge", default=False, action="store_true",
            help="Merge adjacent gaps (to conform to AGP specification)")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    merge = opts.merge
    agpfile, = args
    
    if merge:
        merged_agpfile = agpfile.replace(".agp", ".merged.agp")
        fw = open(merged_agpfile, "w")

    agp = AGP(agpfile)
    size_distribution = defaultdict(int)
    data = [] # store merged AGPLine's

    for is_gap, alines in groupby(agp, key=lambda x: (x.object, x.is_gap)):
        alines = list(alines)
        is_gap = is_gap[1]
        if is_gap:
            gap_size = sum(x.gap_length for x in alines)
            gap_types = set(x.gap_type for x in alines)
            for gtype in ("centromere", "telomere"):
                if gtype in gap_types:
                    gap_size = gtype

            size_distribution[gap_size] += 1
            b = deepcopy(alines[0])
            b.object_beg = min(x.object_beg for x in alines)
            b.object_end = max(x.object_end for x in alines)
            b.gap_length = sum(x.gap_length for x in alines)
        
            assert b.gap_length == b.object_end - b.object_beg + 1
            
            gtypes = [x.gap_type for x in alines]
            for gtype in ("centromere", "telomere", "contig", "clone", "fragment"):
                if gtype in gtypes:
                    b.gap_type = gtype
                    break

            linkages = [x.linkage for x in alines]
            for linkage in ("no", "yes"):
                if linkage in linkages:
                    b.linkage = linkage
                    break

            alines = [b]

        data.extend(alines)

    for gap_size, counts in sorted(size_distribution.items()):
        print >> sys.stderr, gap_size, counts

    if merge:
        AGP.print_header(fw)
        for ob, bb in groupby(data, lambda x: x.object):
            for i, b in enumerate(bb):
                b.part_number = i + 1
                print >> fw, b
    

def build(args):
    """
    %prog build agpfile componentfasta targetfasta
    
    build targetfasta based on info from agpfile
    """
    p = OptionParser(build.__doc__)
    p.add_option("--novalidate", dest="novalidate", default=False,
            action="store_true", 
            help="don't validate the agpfile [default: %default]")

    set_debug(p, args)
    opts, args = p.parse_args(args)

    validate = not (opts.novalidate)

    try:
        agpfile, componentfasta, targetfasta  = args
    except Exception, e:
        sys.exit(p.print_help())

    agp = AGP(agpfile, validate=validate)
    agp.build_all(componentfasta=componentfasta, targetfasta=targetfasta)


def validate(args):
    """
    %prog validate agpfile componentfasta targetfasta
    
    validate consistency between agpfile and targetfasta
    """
    p = OptionParser(validate.__doc__)
    
    set_debug(p, args)
    opts, args = p.parse_args(args)

    try:
        agpfile, componentfasta, targetfasta  = args
    except Exception, e:
        sys.exit(p.print_help())

    agp = AGP(agpfile)
    build = Fasta(targetfasta)
    bacs = Fasta(componentfasta, index=False)
    
    # go through this line by line
    for aline in agp:
        try:
            build_seq = build.sequence(dict(chr=aline.object,
                start=aline.object_beg, stop=aline.object_end))

            if aline.is_gap:
                assert build_seq.upper()==aline.gap_length * 'N', \
                    "gap mismatch: %s" % aline
            else:
                bac_seq = bacs.sequence(dict(chr=aline.component_id,
                    start=aline.component_beg, stop=aline.component_end,
                    strand=aline.orientation))

                assert build_seq.upper()==bac_seq.upper(), \
                        "sequence mismatch: %s" % aline

            logging.debug("%s:%d-%d verified" % (aline.object, aline.object_beg,
                aline.object_end))

        except Exception as e:
            logging.error(str(e))


if __name__ == '__main__':
    main()
