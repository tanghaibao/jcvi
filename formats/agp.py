#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Genbank AGP file format, see spec here
http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp
"""

import os
import re
import sys
import shutil
import logging

from copy import deepcopy
from optparse import OptionParser
from collections import defaultdict
from itertools import groupby

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from jcvi.formats.base import LineFile, must_open
from jcvi.formats.fasta import Fasta
from jcvi.formats.bed import Bed, BedLine
from jcvi.assembly.base import calculate_A50
from jcvi.utils.range import range_intersect
from jcvi.utils.iter import pairwise, flatten
from jcvi.apps.base import ActionDispatcher, set_outfile


Valid_component_type = list("ADFGNOPUW")
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
        self.is_gap = (self.component_type in ('N', 'U'))

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
        return "\t".join((self.object, str(self.object_beg - 1),
                str(self.object_end), gid,
                self.component_type, self.orientation))

    @property
    def isCloneGap(self):
        return (self.is_gap and self.gap_type != "fragment")

    def validate(self):
        assert self.component_type in Valid_component_type, \
                "component_type has to be one of %s" % Valid_component_type
        assert self.object_beg <= self.object_end, \
                "object_beg needs to be <= object_end"

        if not self.is_gap:
            assert self.component_beg <= self.component_end, \
                    "component_begin must be <= component_end"
            assert self.object_span == self.component_span, \
                    "object_span (%d) must be same as component_span (%d)" %\
                    (self.object_span, self.component_span)
        else:
            assert self.gap_length >= 1, \
                    "gap_length must be >= 1"
            assert self.object_span == self.gap_length, \
                    "object span (%d) must be same as gap_length (%d)" % \
                    (self.object_span, self.gap_length)


class AGP (LineFile):

    def __init__(self, filename, validate=True):
        super(AGP, self).__init__(filename)

        fp = open(filename)
        for row in fp:
            if row[0] == '#':
                continue
            self.append(AGPLine(row, validate=validate))

        self.validate = validate
        if validate:
            self.sort(key=lambda x: (x.object, x.object_beg))
            self.validate_all()

    @property
    def order(self):
        """
        Returns a dict with component_id => (i, agpline)
        """
        d = {}
        for (i, x) in enumerate(self):
            if x.is_gap:
                continue
            xid = x.component_id
            d[xid] = (i, x)

            xid = xid.rsplit(".", 1)[0]  # Remove Genbank version
            d[xid] = (i, x)

        return d

    def getAdjacentClone(self, i, south=True):
        """
        Returns the adjacent clone name.
        """
        rr = self[i + 1:] if south else self[i - 1::-1]
        a = self[i]
        for x in rr:
            if x.object != a.object:
                break
            if x.is_gap:
                if x.isCloneGap:
                    return x
                else:
                    continue
            else:
                return x
        return None

    def getNorthSouthClone(self, i):
        """
        Returns the adjacent clone name from both sides.
        """
        north = self.getAdjacentClone(i, south=False)
        south = self.getAdjacentClone(i)
        return north, south

    @classmethod
    def print_header(cls, fw=sys.stdout, organism="Medicago truncatula",
            taxid=3880, source="J. Craig Venter Institute", comment=None):
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

    def report_stats(self, object, bacs, components, scaffold_sizes):
        from jcvi.utils.cbook import human_size

        nbacs = len(bacs)
        nscaffolds = len(scaffold_sizes)
        a50, l50, n50 = calculate_A50(scaffold_sizes)

        print "\t".join(str(x) for x in (object, nbacs,
            components, nscaffolds, n50,
            human_size(l50, precision=2, target="Mb")))

    def summary_one(self, object, lines):
        bacs = set()
        components = 0
        scaffold_sizes = []
        _scaffold_key = lambda x: x.is_gap and \
                x.linkage == "no"

        for is_gap, scaffold in groupby(lines, key=_scaffold_key):
            if is_gap:
                continue

            scaffold = list(scaffold)
            scaffold_size = 0
            for b in scaffold:
                if b.is_gap:
                    scaffold_size += b.gap_length
                else:
                    bacs.add(b.component_id)
                    components += 1
                    scaffold_size += b.component_span

            scaffold_sizes.append(scaffold_size)

        self.report_stats(object, bacs, components, scaffold_sizes)

        return bacs, components, scaffold_sizes

    def summary_all(self):

        all_bacs = set()
        all_scaffold_sizes = []
        all_components = 0
        for ob, lines_with_same_ob in groupby(self, key=lambda x: x.object):
            lines = list(lines_with_same_ob)
            bacs, components, scaffold_sizes = self.summary_one(ob, lines)
            all_components += components
            all_bacs |= bacs
            all_scaffold_sizes.extend(scaffold_sizes)

        self.report_stats("Total", all_bacs, all_components, all_scaffold_sizes)

    def validate_one(self, object, lines):
        object_beg = lines[0].object_beg
        assert object_beg == 1, \
                "object %s must start at 1 (instead of %d)" % \
                (object, object_beg)

        for a, b in pairwise(lines):
            assert b.object_beg - a.object_end == 1, \
                    "lines not continuous coords between:\n%s\n%s" % \
                    (a, b)

    def validate_all(self):
        for ob, lines_with_same_ob in groupby(self, key=lambda x: x.object):
            lines = list(lines_with_same_ob)
            self.validate_one(ob, lines)

    def build_one(self, object, lines, fasta, fw, newagp=None):
        """
        Construct molecule using component fasta sequence
        """
        components = []

        total_bp = 0
        for line in lines:

            if line.is_gap:
                seq = 'N' * line.gap_length
                if newagp:
                    print >> newagp, line
            else:
                seq = fasta.sequence(dict(chr=line.component_id,
                        start=line.component_beg,
                        stop=line.component_end,
                        strand=line.orientation))
                # Check for dangling N's
                if newagp:
                    trimNs(seq, line, newagp)

            components.append(seq)
            total_bp += len(seq)

            if self.validate:
                assert total_bp == line.object_end, \
                        "cumulative base pairs (%d) does not match (%d)" % \
                        (total_bp, line.object_end)

        if not newagp:
            rec = SeqRecord(Seq(''.join(components)), id=object, description="")
            SeqIO.write([rec], fw, "fasta")
            if len(rec) > 1000000:
                logging.debug("Write object %s to `%s`" % (object, fw.name))

    def build_all(self, componentfasta, targetfasta, newagp=None):
        f = Fasta(componentfasta, index=False)
        fw = open(targetfasta, "w")

        for ob, lines_with_same_ob in groupby(self, key=lambda x: x.object):

            lines = list(lines_with_same_ob)
            self.build_one(ob, lines, f, fw, newagp=newagp)


class TPFLine (object):

    def __init__(self, line):
        args = line.split()
        self.component_id = args[0]
        self.object = args[1]
        if self.is_gap:
            self.gap_type = self.component_id

    def __str__(self):
        return "\t".join((self.component_id, self.object_id))

    @property
    def is_gap(self):
        return self.component_id in Valid_gap_type

    @property
    def isCloneGap(self):
        return (self.is_gap and self.gap_type != "fragment")


class TPF (LineFile):

    def __init__(self, filename):
        super(TPF, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            if row[0] == '#':
                continue
            self.append(TPFLine(row))

    def getAdjacentClone(self, i, south=True):
        """
        Returns adjacent clone name, either the line before or after the current
        line.
        """
        rr = self[i + 1:] if south else self[i - 1::-1]
        a = self[i]
        for x in rr:
            if x.object != a.object:
                break
            return x
        return None

    def getNorthSouthClone(self, i):
        """
        Returns adjacent clone name on both sides.
        """
        north = self.getAdjacentClone(i, south=False)
        south = self.getAdjacentClone(i)
        return north, south


def trimNs(seq, line, newagp):
    """
    Test if the sequences contain dangling N's on both sides. This component
    needs to be adjusted to the 'actual' sequence range.
    """
    start, end = line.component_beg, line.component_end
    size = end - start + 1
    leftNs, rightNs = 0, 0
    for s in seq:
        if s in 'nN':
            leftNs += 1
        else:
            break
    for s in seq[::-1]:
        if s in 'nN':
            rightNs += 1
        else:
            break

    if line.orientation == '-':
        trimstart = start + rightNs
        trimend = end - leftNs
    else:
        trimstart = start + leftNs
        trimend = end - rightNs

    trimrange = (trimstart, trimend)
    oldrange = (start, end)

    if trimrange != oldrange:
        logging.error("Range trimmed of N's: {0} => {1}".format(oldrange,
            trimrange))

        if leftNs:
            print >> newagp, "\t".join(str(x) for x in (line.object, 0, 0, 0,
                    'N', leftNs, "fragment", "yes", ""))
        if trimend > trimstart:
            print >> newagp, "\t".join(str(x) for x in (line.object, 0, 0, 0,
                    line.component_type, line.component_id, trimstart, trimend,
                    line.orientation))
        if rightNs and rightNs != size:
            print >> newagp, "\t".join(str(x) for x in (line.object, 0, 0, 0,
                    'N', rightNs, "fragment", "yes", ""))
    else:
        print >> newagp, line


def main():

    actions = (
        ('summary', 'print out a table of scaffold statistics'),
        ('stats', 'print out a report for length of gaps and components'),
        ('phase', 'given genbank file, get the phase for the HTG BAC record'),
        ('bed', 'print out the tiling paths in bed format'),
        ('gaps', 'print out the distribution of gap sizes'),
        ('tpf', 'print out a list of accessions, aka Tiling Path File'),
        ('chr0', 'build AGP file for unplaced sequences'),
        ('mask', 'mask given ranges in components to gaps'),
        ('liftover', 'given ranges in components, get chromosome ranges'),
        ('reindex', 'assume accurate component order, reindex coordinates'),
        ('tidy', 'run trim=>reindex=>merge sequentially'),
        ('build', 'given agp file and component fasta file, build ' +\
                 'pseudomolecule fasta'),
        ('validate', 'given agp file, component and pseudomolecule fasta, ' +\
                     'validate if the build is correct')
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def stats(args):
    """
    %prog stats agpfile

    Print out a report for length of gaps and components.
    """
    p = OptionParser(stats.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    agpfile, = args

    agp = AGP(agpfile)
    gap_lengths = []
    component_lengths = []
    for a in agp:
        span = a.object_span
        if a.is_gap:
            label = a.gap_type
            gap_lengths.append((span, label))
        else:
            label = "{0}:{1}-{2}".format(a.component_id, a.component_beg, \
                   a.component_end)
            component_lengths.append((span, label))
            if span < 50:
                logging.error("component span too small ({0}):\n{1}".\
                    format(span, a))

    table = dict()
    for label, lengths in zip(("Gaps", "Components"),
            (gap_lengths, component_lengths)):

        table[(label, "Min")] = "{0} ({1})".format(*min(lengths))
        table[(label, "Max")] = "{0} ({1})".format(*max(lengths))
        table[(label, "Sum")] = sum(x[0] for x in lengths)

    from jcvi.utils.table import tabulate

    table = tabulate(table)
    print >> sys.stderr, table


def mask(args):
    """
    %prog mask agpfile bedfile

    Mask given ranges in componets to gaps.
    """
    p = OptionParser(mask.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    agpfile, bedfile = args
    agp = AGP(agpfile)
    bed = Bed(bedfile)
    simple_agp = agp.order
    # agp lines to replace original ones, keyed by the component
    agp_fixes = defaultdict(list)

    newagpfile = agpfile.replace(".agp", ".masked.agp")
    logfile = bedfile.replace(".bed", ".masklog")
    fw = open(newagpfile, "w")
    fwlog = open(logfile, "w")

    for component, intervals in bed.sub_beds():
        print >> fwlog, "\n".join(str(x) for x in intervals)
        i, a = simple_agp[component]
        object = a.object
        component_span = a.component_span
        orientation = a.orientation
        print >> fwlog, a

        assert a.component_beg, a.component_end
        arange = a.component_beg, a.component_end

        # Make sure `ivs` contain DISJOINT ranges, and located within `arange`
        ivs = []
        for i in intervals:
            iv = range_intersect(arange, (i.start, i.end))
            if iv is not None:
                ivs.append(iv)

        # Sort the ends of `ivs` as well as the arange
        arange = a.component_beg - 1, a.component_end + 1
        endpoints = sorted(flatten(ivs + [arange]))
        # reverse if component on negative strand
        if orientation == '-':
            endpoints.reverse()

        sum_of_spans = 0
        # assign complements as sequence components
        for i, (a, b) in enumerate(pairwise(endpoints)):
            if orientation == '-':
                a, b = b, a
            if orientation not in ('+', '-'):
                orientation = '+'

            aline = "\t".join(str(x) for x in (object, 0, 0, 0))
            if i % 2 == 0:
                cspan = b - a - 1
                aline = "\t".join(str(x) for x in (aline,
                        'D', component, a + 1, b - 1, orientation))
            else:
                cspan = b - a + 1
                aline = "\t".join(str(x) for x in (aline,
                        "N", cspan, "fragment", "yes"))
            if cspan <= 0:
                continue

            sum_of_spans += cspan
            agp_fixes[component].append(aline)
            print >> fwlog, aline

        assert component_span == sum_of_spans
        print >> fwlog

    # Finally write the masked agp
    for a in agp:
        if not a.is_gap and a.component_id in agp_fixes:
            print >> fw, "\n".join(agp_fixes[a.component_id])
        else:
            print >> fw, a


def liftover(args):
    """
    %prog liftover agpfile bedfile

    Given coordinates in components, convert to the coordinates in chromosomes.
    """
    p = OptionParser(liftover.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    agpfile, bedfile = args
    agp = AGP(agpfile).order
    bed = Bed(bedfile)
    newbed = Bed()
    for b in bed:
        component = b.seqid
        i, a = agp[component]

        assert a.component_beg < a.component_end
        arange = a.component_beg, a.component_end
        assert b.start < b.end
        brange = b.start, b.end

        st = range_intersect(arange, brange)
        if not st:
            continue
        start, end = st
        assert start <= end

        if a.orientation == '-':
            d = a.object_end + a.component_beg
            s, t = d - end, d - start
        else:
            d = a.object_beg - a.component_beg
            s, t = d + start, d + end

        name = "{0}_{1}".format(component, b.accn.replace(" ", "_"))
        bline = "\t".join(str(x) for x in (a.object, s - 1, t, name))
        newbed.append(BedLine(bline))

    newbed.sort(key=lambda x: (x.seqid, x.start))
    newbed.print_to_file()


def reindex(args):
    """
    %prog agpfile

    assume the component line order is correct, modify coordinates, this is
    necessary mostly due to manual edits (insert/delete) that disrupts
    the target coordinates.
    """
    p = OptionParser(reindex.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    agpfile, = args
    agp = AGP(agpfile, validate=False)
    newagpfile = agpfile.replace(".agp", ".reindexed.agp")

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
    logging.error("File `{0}` written and verified.".format(newagpfile))


def summary(args):
    """
    %prog summary agpfile

    print a table of scaffold statistics, number of BACs, no of scaffolds,
    scaffold N50, scaffold L50, actual sequence, PSMOL NNNs, PSMOL-length, % of
    PSMOL sequenced.
    """
    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    agpfile, = args
    header = "Chromosome|# of BACs|# of Components|"
    header += "# of Scaffolds|Scaff N50|Scaff L50"
    header = header.replace('|', '\t')

    agp = AGP(agpfile)
    print header
    agp.summary_all()


chr_pat = re.compile("chromosome (\d)", re.I)
clone_pat = re.compile("clone ([^, ]*\d)[ ,]", re.I)


def get_clone(rec):
    """
    >>> get_clone("Medicago truncatula chromosome 2 clone mth2-48e18")
    ('2', 'mth2-48e18')
    """
    s = rec.description
    chr = re.search(chr_pat, s)
    clone = re.search(clone_pat, s)
    chr = chr.group(1) if chr else ""
    clone = clone.group(1) if clone else ""

    return chr, clone


def get_phase(rec):
    keywords = rec.annotations["keywords"]
    description = rec.description.upper()

    if "HTGS_PHASE1" in keywords:
        phase = 1
    elif "HTGS_PHASE2" in keywords:
        phase = 2
    elif len(keywords) == 1 and "HTG" in keywords:
        phase = 3
    elif "PLN" in keywords:  # EMBL BACs
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
        phase = 3

    return phase, keywords


def phase(args):
    """
    %prog phase genbankfiles

    Input has to be gb file. Search the `KEYWORDS` section to look for PHASE.
    Also look for "chromosome" and "clone" in the definition line.
    """
    p = OptionParser(phase.__doc__)
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fw = must_open(opts.outfile, "w")
    for gbfile in args:
        for rec in SeqIO.parse(gbfile, "gb"):
            bac_phase, keywords = get_phase(rec)
            chr, clone = get_clone(rec)
            keyword_field = ";".join(keywords)
            print >> fw, "\t".join((rec.id, str(bac_phase), keyword_field,
                    chr, clone))


"""
phase 0 - (P)refinish; phase 1,2 - (D)raft;
phase 3 - (F)inished; 4 - (O)thers
"""
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
        comment="{} components with unplaced chr locations".format(len(f)))

    gap_length = opts.gapsize
    object_beg = 1

    if gap_length:
        object = "chr0"
        gap_type = "clone"
        linkage = "no"

        part_number = 0
        for component_id, size in f.itersizes_ordered():
            if part_number > 0:  # print gap except for the first one
                object_end = object_beg + gap_length - 1
                part_number += 1
                print >> fw, "\t".join(str(x) for x in \
                        (object, object_beg, object_end, part_number,
                         'N', gap_length, gap_type, linkage, ""))

                object_beg += gap_length

            object_end = object_beg + size - 1
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


def tpf(args):
    """
    %prog tpf agpfile

    Print out a list of ids, one per line. Also known as the Tiling Path.

    AC225490.9  chr6

    Can optionally output scaffold gaps.
    """
    p = OptionParser(tpf.__doc__)
    p.add_option("--noversion", default=False, action="store_true",
            help="Remove trailing accession versions [default: %default]")
    p.add_option("--gaps", default=False, action="store_true",
            help="Include gaps in the output [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    agpfile, = args
    agp = AGP(agpfile)
    for a in agp:
        object = a.object
        if a.is_gap:
            if opts.gaps and a.isCloneGap:
                print "\t".join((a.gap_type, object))
            continue

        component_id = a.component_id

        if opts.noversion:
            component_id = component_id.rsplit(".", 1)[0]

        print "\t".join((component_id, object))


def bed(args):
    """
    %prog bed agpfile

    print out the tiling paths in bed format
    """
    p = OptionParser(bed.__doc__)
    p.add_option("--gaps", dest="gaps", default=False, action="store_true",
            help="only print bed lines for gaps")
    p.add_option("--nogaps", dest="nogaps", default=False, action="store_true",
            help="do not print bed lines for gaps")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    agpfile, = args
    agp = AGP(agpfile)
    for a in agp:
        if opts.nogaps and a.is_gap:
            continue
        if opts.gaps and not a.is_gap:
            continue
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
    data = []  # store merged AGPLine's

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
            priorities = ("centromere", "telomere", "contig", \
                    "clone", "fragment")
            for gtype in priorities:
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


def tidy(args):
    """
    %prog tidy agpfile componentfasta

    Given an agp file, run through the following steps:
    o Trim components with dangling N's
    o Reindex the agp
    o Merge adjacent gaps

    Final output is in `.tidy.agp`.
    """
    p = OptionParser(tidy.__doc__)

    if len(args) != 2:
        sys.exit(p.print_help())

    agpfile, componentfasta = args
    originalagpfile = agpfile
    tmpfasta = "tmp.fasta"
    build([agpfile, componentfasta, tmpfasta, "--newagp", "--novalidate"])
    os.remove(tmpfasta)

    agpfile = agpfile.replace(".agp", ".trimmed.agp")
    reindex([agpfile])
    os.remove(agpfile)

    agpfile = agpfile.replace(".agp", ".reindexed.agp")
    gaps([agpfile, "--merge"])
    os.remove(agpfile)

    agpfile = agpfile.replace(".agp", ".merged.agp")
    tidyagpfile = originalagpfile.replace(".agp", ".tidy.agp")
    shutil.move(agpfile, tidyagpfile)

    logging.debug("File written to `{0}`.".format(tidyagpfile))


def build(args):
    """
    %prog build agpfile componentfasta targetfasta

    Build targetfasta based on info from agpfile
    """
    p = OptionParser(build.__doc__)
    p.add_option("--newagp", dest="newagp", default=False, action="store_true",
            help="check components to trim dangling N's [default: %default]")
    p.add_option("--novalidate", dest="novalidate", default=False,
            action="store_true",
            help="don't validate the agpfile [default: %default]")

    opts, args = p.parse_args(args)

    validate = not (opts.novalidate)

    if len(args) != 3:
        sys.exit(p.print_help())

    agpfile, componentfasta, targetfasta = args
    if opts.newagp:
        assert agpfile.endswith(".agp")
        newagpfile = agpfile.replace(".agp", ".trimmed.agp")
        newagp = open(newagpfile, "w")
    else:
        newagp = None

    agp = AGP(agpfile, validate=validate)
    agp.build_all(componentfasta=componentfasta, targetfasta=targetfasta,
            newagp=newagp)


def validate(args):
    """
    %prog validate agpfile componentfasta targetfasta

    validate consistency between agpfile and targetfasta
    """
    p = OptionParser(validate.__doc__)

    opts, args = p.parse_args(args)

    try:
        agpfile, componentfasta, targetfasta = args
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
                assert build_seq.upper() == aline.gap_length * 'N', \
                    "gap mismatch: %s" % aline
            else:
                bac_seq = bacs.sequence(dict(chr=aline.component_id,
                    start=aline.component_beg, stop=aline.component_end,
                    strand=aline.orientation))

                assert build_seq.upper() == bac_seq.upper(), \
                        "sequence mismatch: %s" % aline

            logging.debug("%s:%d-%d verified" % (aline.object,
                aline.object_beg, aline.object_end))

        except Exception as e:
            logging.error(e)


if __name__ == '__main__':
    main()
