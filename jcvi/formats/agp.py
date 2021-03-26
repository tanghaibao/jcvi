#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Genbank AGP file format, see spec here
http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp
"""
from __future__ import print_function

import os
import re
import sys
import shutil
import logging

from copy import deepcopy
from collections import defaultdict
from itertools import groupby, zip_longest
from more_itertools import pairwise, flatten

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from jcvi.formats.base import LineFile, must_open
from jcvi.formats.fasta import Fasta
from jcvi.formats.bed import Bed
from jcvi.assembly.base import calculate_A50
from jcvi.utils.range import range_intersect
from jcvi.apps.base import OptionParser, OptionGroup, ActionDispatcher, need_update


Valid_component_type = list("ADFGNOPUW")

Valid_gap_type = (
    "fragment",
    "clone",  # in v1.1, obsolete in v2.0
    "contig",
    "centromere",
    "short_arm",  # in both versions
    "heterochromatin",
    "telomere",
    "repeat",  # in both versions
    "scaffold",
)  # new in v2.0

Valid_orientation = ("+", "-", "0", "?", "na")

Valid_evidence = (
    "",
    "na",
    "paired-ends",
    "align_genus",
    "align_xgenus",
    "align_trnscpt",
    "within_clone",
    "clone_contig",
    "map",
    "strobe",
    "unspecified",
)

component_RGB = {"O": "0,100,0", "F": "0,100,0", "D": "50,205,50", "N": "255,255,255"}

"""
phase 0 - (P)refinish; phase 1,2 - (D)raft;
phase 3 - (F)inished; 4 - (O)thers
"""
Phases = "PDDFO"


class AGPLine(object):
    def __init__(self, row, validate=True):

        atoms = row.split("\t")
        atoms[-1] = atoms[-1].strip()
        self.object = atoms[0]
        self.object_beg = int(atoms[1])
        self.object_end = int(atoms[2])
        self.object_span = self.object_end - self.object_beg + 1
        self.part_number = atoms[3]
        self.component_type = atoms[4]
        self.is_gap = self.component_type in ("N", "U")

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
            self.linkage_evidence = []
            if len(atoms) > 8:
                linkage_evidence = atoms[8].strip()
                if linkage_evidence:
                    self.linkage_evidence = linkage_evidence.split(";")
            self.orientation = "na"
            self.component_id = "{0}.gap{1:03d}".format(
                self.gap_type, int(self.part_number)
            )

        if validate:
            try:
                self.validate()
            except AssertionError as b:
                logging.error("%s\nerror when validating this line:\n%s" % (b, row))

        self.sign = {"+": 1, "-": -1, "?": 0}.get(self.orientation)

    def __str__(self):

        fields = [
            self.object,
            self.object_beg,
            self.object_end,
            self.part_number,
            self.component_type,
        ]

        if not self.is_gap:
            fields += [
                self.component_id,
                self.component_beg,
                self.component_end,
                self.orientation,
            ]
        else:
            fields += [
                self.gap_length,
                self.gap_type,
                self.linkage,
                ";".join(self.linkage_evidence),
            ]

        return "\t".join(str(x) for x in fields)

    __repr__ = __str__

    @property
    def bedline(self):
        # bed formatted line
        gid = self.component_id if not self.is_gap else self.gap_type
        return "\t".join(
            (
                self.object,
                str(self.object_beg - 1),
                str(self.object_end),
                gid,
                self.component_type,
                self.orientation,
            )
        )

    @property
    def bedextra(self):
        # extra lines for bed12
        return "\t".join(
            str(x)
            for x in (
                self.object_beg - 1,
                self.object_end,
                component_RGB[self.component_type],
                1,
                str(self.object_end - self.object_beg + 1) + ",",
                "0,",
            )
        )

    @property
    def bed12line(self):
        # bed12 formatted line
        return self.bedline + "\t" + self.bedextra

    def gffline(self, gff_source="MGSC", gff_feat_type="golden_path_fragment"):
        # gff3 formatted line
        gff_feat_id = "".join(
            str(x) for x in (self.object, ".", "{0:03d}".format(int(self.part_number)))
        )
        attributes = ";".join(
            (
                "ID=" + gff_feat_id,
                "Name=" + self.component_id,
                "phase=" + self.component_type,
            )
        )
        gff_feat_type = "gap" if self.component_type in ["N", "U"] else gff_feat_type
        orientation = "." if self.orientation == "na" else self.orientation

        return "\t".join(
            str(x)
            for x in (
                self.object,
                gff_source,
                gff_feat_type,
                str(self.object_beg),
                str(self.object_end),
                ".",
                orientation,
                ".",
                attributes,
            )
        )

    @property
    def isCloneGap(self):
        return self.is_gap and self.gap_type != "fragment"

    def validate(self):
        assert (
            self.orientation in Valid_orientation
        ), "orientation must be one of {0}".format("|".join(Valid_orientation))
        assert (
            self.component_type in Valid_component_type
        ), "component_type must be one of {0}".format("|".join(Valid_component_type))
        assert (
            self.object_beg <= self.object_end
        ), "object_beg needs to be <= object_end"

        if not self.is_gap:
            assert (
                self.component_beg <= self.component_end
            ), "component_begin must be <= component_end"
            assert (
                self.object_span == self.component_span
            ), "object_span (%d) must be same as component_span (%d)" % (
                self.object_span,
                self.component_span,
            )
        else:
            assert self.gap_length >= 1, "gap_length must be >= 1"
            assert (
                self.object_span == self.gap_length
            ), "object span (%d) must be same as gap_length (%d)" % (
                self.object_span,
                self.gap_length,
            )
            assert all(
                x in Valid_evidence for x in self.linkage_evidence
            ), "linkage_evidence must be one of {0}, you have {1}".format(
                "|".join(Valid_evidence), self.linkage_evidence
            )

            if self.linkage == "no":
                assert not self.linkage_evidence or self.linkage_evidence[0] in (
                    "",
                    "na",
                ), "linkage no is incompatible with evidence {0}".format(
                    self.linkage_evidence
                )

    @classmethod
    def agpline(cls, tuple):
        return AGPLine("\t".join(str(x) for x in tuple), validate=False)

    @classmethod
    def cline(cls, object, cid, sizes, o):
        line = [object, 0, 0, 0]
        line += ["W", cid, 1, sizes[cid], o]
        return AGPLine.agpline(line)

    @classmethod
    def gline(cls, object, gap, unknown=100):
        line = [object, 0, 0, 0]
        gtype = "N"
        if gap < unknown:
            gtype = "U"
            gap = unknown  # Reset it to 100
        line += [gtype, gap, "scaffold", "yes", "paired-ends"]
        return AGPLine.agpline(line)


class AGP(LineFile):
    def __init__(self, filename, nogaps=False, validate=True, sorted=True):
        super(AGP, self).__init__(filename)

        fp = must_open(filename)
        self.header = []
        for row in fp:
            if row[0] == "#":
                self.header.append(row.strip())
                continue
            if row.strip() == "":
                continue
            a = AGPLine(row, validate=validate)
            if nogaps and a.is_gap:
                continue
            self.append(a)

        self.validate = validate
        if validate:
            if not sorted:
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
            if xid not in d:
                d[xid] = (i, x)

        return d

    def getAdjacentClone(self, i, south=True):
        """
        Returns the adjacent clone name.
        """
        rr = range(i + 1, len(self)) if south else range(i - 1, -1, -1)
        a = self[i]
        for ix in rr:
            x = self[ix]
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

    def transfer_header(self, fw=sys.stdout):
        """
        transfer_header() copies header to a new file.
        print_header() creates a new header.
        """
        print("\n".join(self.header), file=fw)

    @classmethod
    def print_header(
        cls, fw=sys.stdout, organism=None, taxid=None, source=None, comment=None
    ):
        # these comments are entirely optional, modeled after maize AGP
        if organism:
            print("# ORGANISM: {0}".format(organism), file=fw)
        if taxid:
            print("# TAX_ID: {0}".format(taxid), file=fw)
        if source:
            print("# GENOME CENTER: {0}".format(source), file=fw)
        if comment:
            print("# COMMENT: {0}".format(comment), file=fw)
        fields = (
            "object object_beg object_end part_number component_type "
            "component_id/gap_length component_beg/gap_type "
            "component_end/linkage orientation/linkage_evidence"
        )
        print("# FIELDS: {0}".format(", ".join(fields.split())), file=fw)

    def rstats(self, object, bacs, components, scaffold_sizes, length):
        from jcvi.utils.cbook import human_size

        nbacs = len(bacs)
        nscaffolds = len(scaffold_sizes)
        a50, l50, n50 = calculate_A50(scaffold_sizes)
        l50 = human_size(l50)
        length = human_size(length)

        return (object, nbacs, components, nscaffolds, n50, l50, length)

    def iter_object(self):
        for ob, lines_with_same_ob in groupby(self, key=lambda x: x.object):
            yield ob, list(lines_with_same_ob)

    def iter_paired_components(self):
        for object, lines in self.iter_object():
            lines = [x for x in lines if not x.is_gap]
            for a, b in pairwise(lines):
                qreverse = a.orientation == "-"
                yield a, b, qreverse

    def print_to_file(self, filename, index=True):
        fw = open(filename, "w")
        for a in self:
            print(a, file=fw)
        fw.close()
        logging.debug("AGP file written to `%s`.", filename)
        if index:
            reindex([filename, "--inplace"])

    def summary_one(self, object, lines):
        bacs = set()
        components = 0
        scaffold_sizes = []
        _scaffold_key = lambda x: x.is_gap and x.linkage == "no"
        length = max(x.object_end for x in lines)

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

        return (
            self.rstats(object, bacs, components, scaffold_sizes, length),
            (bacs, components, scaffold_sizes, length),
        )

    def summary_all(self):

        all_bacs = set()
        all_scaffold_sizes = []
        all_components = 0
        all_length = 0
        for ob, lines in self.iter_object():
            s, bstats = self.summary_one(ob, lines)
            yield s

            bacs, components, scaffold_sizes, length = bstats
            all_components += components
            all_bacs |= bacs
            all_scaffold_sizes.extend(scaffold_sizes)
            all_length += length

        yield self.rstats(
            "Total", all_bacs, all_components, all_scaffold_sizes, all_length
        )

    def validate_one(self, object, lines):
        object_beg = lines[0].object_beg
        assert object_beg == 1, "object %s must start at 1 (instead of %d)" % (
            object,
            object_beg,
        )

        for a, b in pairwise(lines):
            assert (
                b.object_beg - a.object_end == 1
            ), "lines not continuous coords between:\n%s\n%s" % (a, b)

    def validate_all(self):
        for ob, lines in self.iter_object():
            self.validate_one(ob, lines)

    def build_one(self, object, lines, fasta, fw, newagp=None):
        """
        Construct molecule using component fasta sequence
        """
        components = []

        total_bp = 0
        for line in lines:

            if line.is_gap:
                seq = "N" * line.gap_length
                if newagp:
                    print(line, file=newagp)
            else:
                seq = fasta.sequence(
                    dict(
                        chr=line.component_id,
                        start=line.component_beg,
                        stop=line.component_end,
                        strand=line.orientation,
                    )
                )
                # Check for dangling N's
                if newagp:
                    trimNs(seq, line, newagp)

            components.append(seq)
            total_bp += len(seq)

            if self.validate:
                assert (
                    total_bp == line.object_end
                ), "cumulative base pairs (%d) does not match (%d)" % (
                    total_bp,
                    line.object_end,
                )

        if not newagp:
            rec = SeqRecord(Seq("".join(components)), id=object, description="")
            SeqIO.write([rec], fw, "fasta")
            if len(rec) > 1000000:
                logging.debug("Write object %s to `%s`" % (object, fw.name))

    def build_all(self, componentfasta, targetfasta, newagp=None):
        f = Fasta(componentfasta, index=False)
        fw = open(targetfasta, "w")

        for ob, lines in self.iter_object():
            self.build_one(ob, lines, f, fw, newagp=newagp)

    @property
    def graph(self):
        from jcvi.algorithms.graph import BiGraph

        g = BiGraph()
        for ob, lines in self.iter_object():
            components = [x for x in lines if not x.is_gap]
            gaps = [x for x in lines if x.is_gap]
            for i, (a, b) in enumerate(pairwise(components)):
                g.add_edge(
                    a.component_id,
                    b.component_id,
                    a.orientation,
                    b.orientation,
                    length=gaps[i].gap_length,
                )
            if len(components) == 1:  # Singleton object
                a = components[0]
                g.add_node(a.component_id)

        return g

    def get_line(self, cid):
        for i, a in enumerate(self):
            if not a.is_gap and a.component_id == cid:
                return i, a
        return None, None

    # Update AGP on the fly
    def delete_line(self, a, verbose=False):
        ai, ax = self.get_line(a)
        if ai is None:
            return

        if verbose:
            msg = "* Delete line:\n{0}".format(ax)
            print(msg, file=sys.stderr)

        del self[ai]

    def delete_lines(self, lines, verbose=False):
        deleted = set()
        for r in lines:
            if r.is_gap:
                continue
            cid = r.component_id
            self.delete_line(cid, verbose=verbose)
            deleted.add(cid)
        return deleted

    def insert_lines(self, a, lines, after=False, delete=False, verbose=False):
        if delete:
            deleted = self.delete_lines(lines, verbose=verbose)

        ai, ax = self.get_line(a)
        if after:
            ai += 1
        for i, x in enumerate(lines):
            self.insert(ai + i, x)
        if verbose:
            tag = "after" if after else "before"
            msg = "* Insert {0} line:\n".format(tag)
            msg += "\n".join([str(ax), "-" * 60]) + "\n"
            msg += "\n".join(str(x) for x in lines)
            print(msg, file=sys.stderr)
        return deleted

    def update_between(self, a, b, lines, delete=True, verbose=False):
        if delete:
            deleted = self.delete_lines(lines, verbose=verbose)

        ai, ax = self.get_line(a)
        bi, bx = self.get_line(b)
        # Update
        self[ai + 1 : bi] = lines
        if verbose:
            msg = "* Update between:\n"
            msg += "\n".join([str(ax), str(bx), "-" * 60]) + "\n"
            msg += "\n".join(str(x) for x in lines)
            print(msg, file=sys.stderr)
        return deleted

    def convert_to_gap(self, a, verbose=False):
        ai, ax = self.get_line(a)
        gline = AGPLine.gline(ax.object, 100)
        self[ai] = gline
        if verbose:
            msg = "* Convert from/to:\n"
            msg += "\n".join([str(ax), str(gline), "-" * 60]) + "\n"
            print(msg, file=sys.stderr)

    def delete_between(self, a, b, verbose=True):
        return self.update_between(a, b, [], verbose=verbose)

    def switch_between(self, a, b, verbose=True):
        ai, ax = self.get_line(a)
        bi, bx = self.get_line(b)
        self[ai] = bx
        self[bi] = ax
        if verbose:
            msg = "* Switch between:\n"
            msg += "\n".join([str(ax), str(bx)])
            print(msg, file=sys.stderr)


class TPFLine(object):
    def __init__(self, line):
        args = line.split()
        self.component_id = args[0]
        self.object = args[1]
        if self.is_gap:
            self.gap_type = self.component_id
        self.orientation = args[2]

    def __str__(self):
        return "\t".join((self.component_id, self.object_id, self.orientation))

    @property
    def is_gap(self):
        return self.component_id in Valid_gap_type

    @property
    def isCloneGap(self):
        return self.is_gap and self.gap_type != "fragment"


class TPF(LineFile):
    def __init__(self, filename):
        super(TPF, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            if row[0] == "#":
                continue
            self.append(TPFLine(row))

    def getAdjacentClone(self, i, south=True):
        """
        Returns adjacent clone name, either the line before or after the current
        line.
        """
        rr = range(i + 1, len(self)) if south else range(i - 1, -1, -1)
        a = self[i]
        for ix in rr:
            x = self[ix]
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


class OOLine(object):
    def __init__(self, id, component_id, component_size, strand):
        self.id = id
        self.component_id = component_id
        self.component_size = component_size
        self.strand = strand


class OO(LineFile):
    def __init__(self, filename=None, ctgsizes=None):
        super(OO, self).__init__(filename)

        if filename is None:
            return

        from jcvi.formats.base import read_block

        fp = open(filename)
        prefix = "contig_"
        self.contigs = set()
        for header, block in read_block(fp, ">"):
            header = header[1:]  # Trim the '>'
            header = header.split()[0]
            for b in block:
                ctg, orientation = b.split()
                if ctg.startswith(prefix):
                    ctg = ctg[len(prefix) :]

                assert orientation in ("BE", "EB")

                strand = "+" if orientation == "BE" else "-"
                ctgsize = ctgsizes[ctg]
                self.add(header, ctg, ctgsize, strand)
                self.contigs.add(ctg)

    def add(self, scaffold, ctg, ctgsize, strand="0"):
        self.append(OOLine(scaffold, ctg, ctgsize, strand))

    def sub_beds(self):
        for scaffold, beds in groupby(self, key=lambda x: x.id):
            yield scaffold, list(beds)

    def write_AGP(
        self, fw=sys.stdout, gapsize=100, phases={}, gaptype="scaffold", evidence=""
    ):

        linkage = "yes"

        for object, beds in self.sub_beds():
            object_beg = 1
            part_number = 0
            for b in beds:
                component_id = b.component_id
                size = b.component_size
                if (
                    part_number > 0 and gapsize > 0
                ):  # Print gap except for the first one
                    object_end = object_beg + gapsize - 1
                    part_number += 1
                    component_type = "U" if gapsize == 100 else "N"
                    print(
                        "\t".join(
                            str(x)
                            for x in (
                                object,
                                object_beg,
                                object_end,
                                part_number,
                                component_type,
                                gapsize,
                                gaptype,
                                linkage,
                                evidence,
                            )
                        ),
                        file=fw,
                    )

                    object_beg += gapsize

                object_end = object_beg + size - 1
                part_number += 1
                strand = "?" if b.strand == "0" else b.strand
                print(
                    "\t".join(
                        str(x)
                        for x in (
                            object,
                            object_beg,
                            object_end,
                            part_number,
                            phases.get(component_id, "W"),
                            component_id,
                            1,
                            size,
                            strand,
                        )
                    ),
                    file=fw,
                )

                object_beg += size


def order_to_agp(
    object, ctgorder, sizes, fwagp, gapsize=100, gaptype="scaffold", evidence=""
):

    o = OO()  # Without a filename
    for scaffold_number, (ctg, strand) in enumerate(ctgorder):
        size = sizes[ctg]
        o.add(object, ctg, size, strand)

    o.write_AGP(fwagp, gapsize=gapsize, gaptype=gaptype, phases={}, evidence=evidence)


def trimNs(seq, line, newagp):
    """
    Test if the sequences contain dangling N's on both sides. This component
    needs to be adjusted to the 'actual' sequence range.
    """
    start, end = line.component_beg, line.component_end
    size = end - start + 1
    leftNs, rightNs = 0, 0
    lid, lo = line.component_id, line.orientation
    for s in seq:
        if s in "nN":
            leftNs += 1
        else:
            break
    for s in seq[::-1]:
        if s in "nN":
            rightNs += 1
        else:
            break

    if lo == "-":
        trimstart = start + rightNs
        trimend = end - leftNs
    else:
        trimstart = start + leftNs
        trimend = end - rightNs

    trimrange = (trimstart, trimend)
    oldrange = (start, end)

    if trimrange != oldrange:
        logging.debug("{0} trimmed of N's: {1} => {2}".format(lid, oldrange, trimrange))

        if leftNs:
            print(
                "\t".join(
                    str(x)
                    for x in (line.object, 0, 0, 0, "N", leftNs, "fragment", "yes", "")
                ),
                file=newagp,
            )
        if trimend > trimstart:
            print(
                "\t".join(
                    str(x)
                    for x in (
                        line.object,
                        0,
                        0,
                        0,
                        line.component_type,
                        lid,
                        trimstart,
                        trimend,
                        lo,
                    )
                ),
                file=newagp,
            )
        if rightNs and rightNs != size:
            print(
                "\t".join(
                    str(x)
                    for x in (line.object, 0, 0, 0, "N", rightNs, "fragment", "yes", "")
                ),
                file=newagp,
            )
    else:
        print(line, file=newagp)


def main():

    actions = (
        ("summary", "print out a table of scaffold statistics"),
        ("stats", "print out a report for length of gaps and components"),
        ("phase", "given genbank file, get the phase for the HTG BAC record"),
        ("bed", "print out the tiling paths in bed/gff3 format"),
        ("frombed", "generate AGP file based on bed file"),
        ("fromcsv", "generate AGP file based on simple csv file"),
        (
            "extendbed",
            "extend the components to fill the component range and output bed/gff3 format file",
        ),
        ("gaps", "print out the distribution of gap sizes"),
        ("tpf", "print out a list of accessions, aka Tiling Path File"),
        ("cut", "cut at the boundaries of given ranges"),
        ("mask", "mask given ranges in components to gaps"),
        ("swap", "swap objects and components"),
        ("format", "reformat AGP file"),
        ("reindex", "assume accurate component order, reindex coordinates"),
        ("tidy", "run trim=>reindex=>merge sequentially"),
        (
            "build",
            "given agp file and component fasta file, build " + "pseudomolecule fasta",
        ),
        (
            "validate",
            "given agp file, component and pseudomolecule fasta, "
            + "validate if the build is correct",
        ),
        ("infer", "infer where the components are in the genome"),
        ("compress", "compress coordinates based on multiple AGP files"),
    )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fromcsv(args):
    """
    %prog fromcsv contigs.fasta map.csv map.agp

    Convert csv which contains list of scaffolds/contigs to AGP file.
    """
    import csv
    from jcvi.formats.sizes import Sizes

    p = OptionParser(fromcsv.__doc__)
    p.add_option("--evidence", default="map", help="Linkage evidence to add in AGP")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    contigsfasta, mapcsv, mapagp = args
    reader = csv.reader(open(mapcsv))
    sizes = Sizes(contigsfasta).mapping
    next(reader)  # Header
    fwagp = must_open(mapagp, "w")
    o = OO()
    for row in reader:
        if len(row) == 2:
            object, ctg = row
            strand = "?"
        elif len(row) == 3:
            object, ctg, strand = row
        size = sizes[ctg]
        o.add(object, ctg, size, strand)

    o.write_AGP(
        fwagp, gapsize=100, gaptype="scaffold", phases={}, evidence=opts.evidence
    )


def compress(args):
    """
    %prog compress a.agp b.agp

    Convert coordinates based on multiple AGP files. Useful to simplify multiple
    liftOvers to compress multiple chain files into a single chain file, in
    upgrading locations of genomic features.

    Example:
    `a.agp` could contain split scaffolds:
    scaffold_0.1    1       600309  1       W       scaffold_0      1 600309  +

    `b.agp` could contain mapping to chromosomes:
    LG05    6435690 7035998 53      W       scaffold_0.1    1       600309  +

    The final AGP we want is:
    LG05    6435690 7035998 53      W       scaffold_0      1       600309  +
    """
    p = OptionParser(compress.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    aagpfile, bagpfile = args
    # First AGP provides the mapping
    store = {}
    agp = AGP(aagpfile)
    for a in agp:
        if a.is_gap:
            continue
        # Ignore '?' in the mapping
        if a.sign == 0:
            a.sign = 1
        store[(a.object, a.object_beg, a.object_end)] = (
            a.component_id,
            a.component_beg,
            a.component_end,
            a.sign,
        )

    # Second AGP forms the backbone
    agp = AGP(bagpfile)
    fw = must_open(opts.outfile, "w")
    print("\n".join(agp.header), file=fw)
    for a in agp:
        if a.is_gap:
            print(a, file=fw)
            continue
        component_id, component_beg, component_end, sign = store[
            (a.component_id, a.component_beg, a.component_end)
        ]

        orientation = {1: "+", -1: "-", 0: "?"}.get(sign * a.sign)
        atoms = (
            a.object,
            a.object_beg,
            a.object_end,
            a.part_number,
            a.component_type,
            component_id,
            component_beg,
            component_end,
            orientation,
        )
        a = AGPLine("\t".join(str(x) for x in atoms))
        print(a, file=fw)


def map_one_scaffold_1way(scaffold_name, scaffold, genome, orientation="+"):
    if orientation == "-":
        scaffold = scaffold.reverse_complement()

    scaffold = str(scaffold)
    for obj_name, obj in genome.iteritems():
        obj_idx = obj.find(scaffold)
        if obj_idx == -1:
            continue
        else:
            return obj_name, obj_idx, orientation
    return -1, -1, orientation  # unmapped scaffolds


def map_one_scaffold(opts):
    scaffold_name, scaffold, genome = opts
    scaffold = scaffold.seq
    obj_name, obj_idx, objo = map_one_scaffold_1way(scaffold_name, scaffold, genome)
    if obj_name == -1:
        obj_name, obj_idx, objo = map_one_scaffold_1way(
            scaffold_name, scaffold, genome, orientation="-"
        )
    if obj_name == -1:
        return ""

    obj_end = obj_idx + len(scaffold)
    return "\t".join(
        str(x) for x in (obj_name, obj_idx, obj_end, scaffold_name, 1000, objo)
    )


def check_seen(r, seen):
    from jcvi.utils.range import range_overlap

    for s in seen:
        if range_overlap(r, s):
            return True
    return False


def infer(args):
    """
    %prog infer scaffolds.fasta genome.fasta

    Infer where the components are in the genome. This function is rarely used,
    but can be useful when distributor does not ship an AGP file.
    """
    from jcvi.apps.grid import WriteJobs
    from jcvi.formats.bed import sort

    p = OptionParser(infer.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    scaffoldsf, genomef = args
    inferbed = "infer-components.bed"
    if need_update((scaffoldsf, genomef), inferbed):
        scaffolds = Fasta(scaffoldsf, lazy=True)
        genome = Fasta(genomef)
        genome = genome.tostring()
        args = [
            (scaffold_name, scaffold, genome)
            for scaffold_name, scaffold in scaffolds.iteritems_ordered()
        ]

        pool = WriteJobs(map_one_scaffold, args, inferbed, cpus=opts.cpus)
        pool.run()

    sort([inferbed, "-i"])
    bed = Bed(inferbed)
    inferagpbed = "infer.bed"
    fw = open(inferagpbed, "w")
    seen = []
    for b in bed:
        r = (b.seqid, b.start, b.end)
        if check_seen(r, seen):
            continue
        print(
            "\t".join(str(x) for x in (b.accn, 0, b.span, b.seqid, b.score, b.strand)),
            file=fw,
        )
        seen.append(r)
    fw.close()

    frombed([inferagpbed])


def format(args):
    """
    %prog format oldagpfile newagpfile

    Reformat AGP file. --switch will replace the ids in the AGP file.
    """
    from jcvi.formats.base import DictFile

    p = OptionParser(format.__doc__)
    p.add_option("--switchcomponent", help="Switch component id based on")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    oldagpfile, newagpfile = args
    switchcomponent = opts.switchcomponent
    if switchcomponent:
        switchcomponent = DictFile(switchcomponent)

    agp = AGP(oldagpfile)
    fw = open(newagpfile, "w")
    nconverts = 0
    for i, a in enumerate(agp):
        if not a.is_gap and a.component_id in switchcomponent:
            oldid = a.component_id
            newid = switchcomponent[a.component_id]
            a.component_id = newid
            logging.debug("Covert {0} to {1} on line {2}".format(oldid, newid, i + 1))
            nconverts += 1
        print(a, file=fw)

    logging.debug("Total converted records: {0}".format(nconverts))


def frombed(args):
    """
    %prog frombed bedfile

    Generate AGP file based on bed file. The bed file must have at least 6
    columns. With the 4-th column indicating the new object.
    """
    p = OptionParser(frombed.__doc__)
    p.add_option(
        "--gapsize",
        default=100,
        type="int",
        help="Insert gaps of size",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    gapsize = opts.gapsize
    agpfile = bedfile.replace(".bed", ".agp")
    fw = open(agpfile, "w")

    bed = Bed(bedfile, sorted=False)
    for object, beds in groupby(bed, key=lambda x: x.accn):
        beds = list(beds)
        for i, b in enumerate(beds):
            if gapsize and i != 0:
                print(
                    "\t".join(
                        str(x)
                        for x in (
                            object,
                            0,
                            0,
                            0,
                            "U",
                            gapsize,
                            "scaffold",
                            "yes",
                            "map",
                        )
                    ),
                    file=fw,
                )

            print(
                "\t".join(
                    str(x)
                    for x in (object, 0, 0, 0, "W", b.seqid, b.start, b.end, b.strand)
                ),
                file=fw,
            )

    fw.close()

    # Reindex
    return reindex([agpfile, "--inplace"])


def swap(args):
    """
    %prog swap agpfile

    Swap objects and components. Will add gap lines. This is often used in
    conjuction with formats.chain.fromagp() to convert between different
    coordinate systems.
    """
    from jcvi.utils.range import range_interleave

    p = OptionParser(swap.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (agpfile,) = args

    agp = AGP(agpfile, nogaps=True, validate=False)
    agp.sort(key=lambda x: (x.component_id, x.component_beg))

    newagpfile = agpfile.rsplit(".", 1)[0] + ".swapped.agp"
    fw = open(newagpfile, "w")
    agp.transfer_header(fw)
    for cid, aa in groupby(agp, key=(lambda x: x.component_id)):
        aa = list(aa)
        aranges = [(x.component_id, x.component_beg, x.component_end) for x in aa]
        gaps = range_interleave(aranges)
        for a, g in zip_longest(aa, gaps):
            a.object, a.component_id = a.component_id, a.object
            a.component_beg = a.object_beg
            a.component_end = a.object_end
            print(a, file=fw)
            if not g:
                continue

            aline = [cid, 0, 0, 0]
            gseq, ga, gb = g
            cspan = gb - ga + 1
            aline += ["N", cspan, "fragment", "yes"]
            print("\t".join(str(x) for x in aline), file=fw)

    fw.close()
    # Reindex
    idxagpfile = reindex([newagpfile, "--inplace"])

    return newagpfile


def stats(args):
    """
    %prog stats agpfile

    Print out a report for length of gaps and components.
    """
    from jcvi.utils.table import tabulate

    p = OptionParser(stats.__doc__)
    p.add_option(
        "--warn",
        default=False,
        action="store_true",
        help="Warnings on small component spans",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (agpfile,) = args

    agp = AGP(agpfile)
    gap_lengths = []
    component_lengths = []
    for a in agp:
        span = a.object_span
        if a.is_gap:
            label = a.gap_type
            gap_lengths.append((span, label))
        else:
            label = "{0}:{1}-{2}".format(
                a.component_id, a.component_beg, a.component_end
            )
            component_lengths.append((span, label))
            if opts.warn and span < 50:
                logging.error("component span too small ({0}):\n{1}".format(span, a))

    table = dict()
    for label, lengths in zip(("Gaps", "Components"), (gap_lengths, component_lengths)):

        if not lengths:
            table[(label, "Min")] = table[(label, "Max")] = table[
                (label, "Sum")
            ] = "n.a."
            continue

        table[(label, "Min")] = "{0} ({1})".format(*min(lengths))
        table[(label, "Max")] = "{0} ({1})".format(*max(lengths))
        table[(label, "Sum")] = sum(x[0] for x in lengths)

    print(tabulate(table), file=sys.stderr)


def cut(args):
    """
    %prog cut agpfile bedfile

    Cut at the boundaries of the ranges in the bedfile.
    """
    p = OptionParser(cut.__doc__)
    p.add_option("--sep", default=".", help="Separator for splits")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    agpfile, bedfile = args
    sep = opts.sep

    agp = AGP(agpfile)
    bed = Bed(bedfile)
    simple_agp = agp.order
    newagpfile = agpfile.replace(".agp", ".cut.agp")
    fw = open(newagpfile, "w")

    agp_fixes = defaultdict(list)
    for component, intervals in bed.sub_beds():
        i, a = simple_agp[component]
        object = a.object
        component_span = a.component_span
        orientation = a.orientation

        assert a.component_beg, a.component_end
        cuts = set()
        for i in intervals:
            start, end = i.start, i.end
            end -= 1

            assert start <= end
            cuts.add(start)
            cuts.add(end)

        cuts.add(0)
        cuts.add(component_span)
        cuts = list(sorted(cuts))

        sum_of_spans = 0
        for i, (a, b) in enumerate(pairwise(cuts)):
            oid = object + "{0}{1}".format(sep, i + 1)
            aline = [oid, 0, 0, 0]
            cspan = b - a
            aline += ["D", component, a + 1, b, orientation]
            sum_of_spans += cspan

            aline = "\t".join(str(x) for x in aline)
            agp_fixes[component].append(aline)

        assert component_span == sum_of_spans

    # Finally write the masked agp
    for a in agp:
        if not a.is_gap and a.component_id in agp_fixes:
            print("\n".join(agp_fixes[a.component_id]), file=fw)
        else:
            print(a, file=fw)

    fw.close()
    # Reindex
    reindex([newagpfile, "--inplace"])

    return newagpfile


def mask(args):
    """
    %prog mask agpfile bedfile

    Mask given ranges in components to gaps. When the bedfile contains a single
    base pair, this position can be a point of split and no base is lost
    (--splitsingle).
    """
    p = OptionParser(mask.__doc__)
    p.add_option(
        "--splitobject",
        default=False,
        action="store_true",
        help="Create new names for object",
    )
    p.add_option(
        "--splitcomponent",
        default=False,
        action="store_true",
        help="Create new names for component",
    )
    p.add_option(
        "--splitsingle",
        default=False,
        action="store_true",
        help="Do not remove base on single point",
    )
    p.add_option(
        "--gaptype",
        default="scaffold",
        help="Masked region has gap type of",
    )
    p.add_option(
        "--noretain",
        default=False,
        action="store_true",
        help="Do not retain old names for non-split objects",
    )
    p.add_option("--sep", default=".", help="Separator for splits")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    agpfile, bedfile = args
    gaptype = opts.gaptype
    splitobject = opts.splitobject
    splitcomponent = opts.splitcomponent
    sep = opts.sep

    assert not (
        splitobject and splitcomponent
    ), "Options --splitobject and --splitcomponent conflict"

    agp = AGP(agpfile)
    bed = Bed(bedfile)
    simple_agp = agp.order
    # agp lines to replace original ones, keyed by the component
    agp_fixes = defaultdict(list)

    newagpfile = agpfile.replace(".agp", ".masked.agp")
    fw = open(newagpfile, "w")

    if splitcomponent:
        componentindex = defaultdict(int)

    for component, intervals in bed.sub_beds():
        i, a = simple_agp[component]
        object = a.object
        orientation = a.orientation

        assert a.component_beg, a.component_end
        arange = a.component_beg, a.component_end

        # Make sure `ivs` contain DISJOINT ranges, and located within `arange`
        ivs = []
        points = set()
        for i in intervals:
            start, end = i.start, i.end
            if opts.splitsingle:
                points.add(start)
            iv = range_intersect(arange, (start, end))
            if iv is not None:
                ivs.append(iv)

        # Sort the ends of `ivs` as well as the arange
        arange = a.component_beg - 1, a.component_end + 1
        endpoints = sorted(flatten(ivs + [arange]))
        # reverse if component on negative strand
        if orientation == "-":
            endpoints.reverse()

        sum_of_spans = 0
        # assign complements as sequence components
        for i, (a, b) in enumerate(pairwise(endpoints)):
            if orientation == "-":
                a, b = b, a
            if orientation not in ("+", "-"):
                orientation = "+"

            oid = object + "{0}{1}".format(sep, i // 2 + 1) if splitobject else object
            aline = [oid, 0, 0, 0]
            if i % 2 == 0:
                cspan = b - a - 1
                if splitcomponent:
                    cid = component + "{0}{1}".format(
                        sep, componentindex[component] + 1
                    )
                    componentindex[component] += 1
                    aline += ["W", cid, 1, cspan, orientation]
                else:
                    end = b if (opts.splitsingle and b in points) else b - 1
                    aline += ["W", component, a + 1, end, orientation]
                is_gap = False
            else:
                cspan = b - a + 1
                aline += ["N", cspan, gaptype, "yes", "paired-ends"]
                is_gap = True
            if cspan <= 0:
                continue

            sum_of_spans += cspan
            aline = "\t".join(str(x) for x in aline)
            if not (splitobject and is_gap):
                agp_fixes[component].append(aline)

    retain = not opts.noretain
    # Finally write the masked agp
    for a in agp:
        if a.is_gap:
            print(a, file=fw)
        elif a.component_id in agp_fixes:
            print("\n".join(agp_fixes[a.component_id]), file=fw)
        else:
            if not retain:
                if splitobject:
                    a.object += sep + "0"
                elif splitcomponent:
                    a.component_id += sep + "0"
            print(a, file=fw)

    fw.close()

    # Reindex
    idxagpfile = reindex([newagpfile, "--inplace"])

    return newagpfile


def reindex(args):
    """
    %prog agpfile

    assume the component line order is correct, modify coordinates, this is
    necessary mostly due to manual edits (insert/delete) that disrupts
    the target coordinates.
    """
    p = OptionParser(reindex.__doc__)
    p.add_option(
        "--nogaps",
        default=False,
        action="store_true",
        help="Remove all gap lines",
    )
    p.add_option(
        "--inplace",
        default=False,
        action="store_true",
        help="Replace input file",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (agpfile,) = args
    inplace = opts.inplace
    agp = AGP(agpfile, validate=False)
    pf = agpfile.rsplit(".", 1)[0]
    newagpfile = pf + ".reindexed.agp"

    fw = open(newagpfile, "w")
    agp.transfer_header(fw)
    for chr, chr_agp in groupby(agp, lambda x: x.object):
        chr_agp = list(chr_agp)
        object_beg = 1
        for i, b in enumerate(chr_agp):
            b.object_beg = object_beg
            b.part_number = i + 1
            if opts.nogaps and b.is_gap:
                continue

            if b.is_gap:
                b.object_end = object_beg + b.gap_length - 1
            else:
                b.object_end = object_beg + b.component_span - 1

            object_beg = b.object_end + 1

            print(str(b), file=fw)

    # Last step: validate the new agpfile
    fw.close()
    agp = AGP(newagpfile, validate=True)

    if inplace:
        shutil.move(newagpfile, agpfile)
        logging.debug("Rename file `{0}` to `{1}`".format(newagpfile, agpfile))
        newagpfile = agpfile

    return newagpfile


def summary(args):
    """
    %prog summary agpfile

    print a table of scaffold statistics, number of BACs, no of scaffolds,
    scaffold N50, scaffold L50, actual sequence, PSMOL NNNs, PSMOL-length, % of
    PSMOL sequenced.
    """
    from jcvi.utils.table import write_csv

    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (agpfile,) = args
    header = (
        "Chromosome #_Distinct #_Components #_Scaffolds "
        "Scaff_N50 Scaff_L50 Length".split()
    )

    agp = AGP(agpfile)
    data = list(agp.summary_all())
    write_csv(header, data, sep=" ")


chr_pat = re.compile(r"chromosome (\d)", re.I)
clone_pat = re.compile(r"clone ([^, ]*\d)[ ,]", re.I)


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
        # logging.error("{0}: {1}".format(rec.name, description))
        phase = 3

    return phase, keywords


def phase(args):
    """
    %prog phase genbankfiles

    Input has to be gb file. Search the `KEYWORDS` section to look for PHASE.
    Also look for "chromosome" and "clone" in the definition line.
    """
    p = OptionParser(phase.__doc__)
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fw = must_open(opts.outfile, "w")
    for gbfile in args:
        for rec in SeqIO.parse(gbfile, "gb"):
            bac_phase, keywords = get_phase(rec)
            chr, clone = get_clone(rec)
            keyword_field = ";".join(keywords)
            print(
                "\t".join((rec.id, str(bac_phase), keyword_field, chr, clone)), file=fw
            )


def tpf(args):
    """
    %prog tpf agpfile

    Print out a list of ids, one per line. Also known as the Tiling Path.

    AC225490.9  chr6

    Can optionally output scaffold gaps.
    """
    p = OptionParser(tpf.__doc__)
    p.add_option(
        "--noversion",
        default=False,
        action="store_true",
        help="Remove trailing accession versions",
    )
    p.add_option(
        "--gaps",
        default=False,
        action="store_true",
        help="Include gaps in the output",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (agpfile,) = args
    agp = AGP(agpfile)
    for a in agp:
        object = a.object
        if a.is_gap:
            if opts.gaps and a.isCloneGap:
                print("\t".join((a.gap_type, object, "na")))
            continue

        component_id = a.component_id
        orientation = a.orientation

        if opts.noversion:
            component_id = component_id.rsplit(".", 1)[0]

        print("\t".join((component_id, object, orientation)))


def bed(args):
    """
    %prog bed agpfile

    print out the tiling paths in bed/gff3 format
    """
    from jcvi.formats.obo import validate_term

    p = OptionParser(bed.__doc__)
    p.add_option(
        "--gaps",
        default=False,
        action="store_true",
        help="Only print bed lines for gaps",
    )
    p.add_option(
        "--nogaps",
        default=False,
        action="store_true",
        help="Do not print bed lines for gaps",
    )
    p.add_option(
        "--bed12",
        default=False,
        action="store_true",
        help="Produce bed12 formatted output",
    )
    p.add_option(
        "--component",
        default=False,
        action="store_true",
        help="Generate bed file for components",
    )
    p.set_outfile()
    g1 = OptionGroup(
        p,
        "GFF specific parameters",
        "Note: If not specified, output will be in `bed` format",
    )
    g1.add_option(
        "--gff",
        default=False,
        action="store_true",
        help="Produce gff3 formatted output. By default, ignores AGP gap lines",
    )
    g1.add_option("--source", default="MGSC", help="Specify a gff3 source")
    g1.add_option(
        "--feature",
        default="golden_path_fragment",
        help="Specify a gff3 feature type",
    )
    p.add_option_group(g1)
    p.set_SO_opts()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    if opts.component:
        opts.nogaps = True

    # If output format is gff3 and 'verifySO' option is invoked, validate the SO term
    if opts.gff and opts.verifySO:
        validate_term(opts.feature, method=opts.verifySO)

    (agpfile,) = args
    agp = AGP(agpfile)
    fw = must_open(opts.outfile, "w")
    if opts.gff:
        print("##gff-version 3", file=fw)

    for a in agp:
        if opts.nogaps and a.is_gap:
            continue
        if opts.gaps and not a.is_gap:
            continue
        if opts.bed12:
            print(a.bed12line, file=fw)
        elif opts.gff:
            print(
                a.gffline(gff_source=opts.source, gff_feat_type=opts.feature), file=fw
            )
        elif opts.component:
            name = "{0}:{1}-{2}".format(
                a.component_id, a.component_beg, a.component_end
            )
            print(
                "\t".join(
                    str(x)
                    for x in (
                        a.component_id,
                        a.component_beg - 1,
                        a.component_end,
                        name,
                        a.component_type,
                        a.orientation,
                    )
                ),
                file=fw,
            )
        else:
            print(a.bedline, file=fw)
    fw.close()

    return fw.name


def extendbed(args):
    """
    %prog extend agpfile componentfasta

    Extend the components to fill the component range. For example, a bed/gff3 file
    that was converted from the agp will contain only the BAC sequence intervals
    that are 'represented' - sometimes leaving the 5` and 3` out (those that
    overlap with adjacent sequences. This script fill up those ranges,
    potentially to make graphics for tiling path.
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(extendbed.__doc__)
    p.add_option(
        "--nogaps",
        default=False,
        action="store_true",
        help="Do not print bed lines for gaps",
    )
    p.add_option(
        "--bed12",
        default=False,
        action="store_true",
        help="Produce bed12 formatted output",
    )
    p.add_option(
        "--gff",
        default=False,
        action="store_true",
        help="Produce gff3 formatted output. By default, ignores " + " AGP gap lines.",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    # If output format is GFF3, ignore AGP gap lines.
    if opts.gff:
        opts.nogaps = True

    agpfile, fastafile = args
    agp = AGP(agpfile)
    fw = must_open(opts.outfile, "w")
    if opts.gff:
        print("##gff-version 3", file=fw)

    ranges = defaultdict(list)
    thickCoords = []  # These are the coordinates before modify ranges
    # Make the first pass to record all the component ranges
    for a in agp:
        thickCoords.append((a.object_beg, a.object_end))
        if a.is_gap:
            continue
        ranges[a.component_id].append(a)

    # Modify the ranges
    sizes = Sizes(fastafile).mapping
    for accn, rr in ranges.items():
        alen = sizes[accn]

        a = rr[0]
        if a.orientation == "+":
            hang = a.component_beg - 1
        else:
            hang = alen - a.component_end
        a.object_beg -= hang

        a = rr[-1]
        if a.orientation == "+":
            hang = alen - a.component_end
        else:
            hang = a.component_beg - 1
        a.object_end += hang

    for a, (ts, te) in zip(agp, thickCoords):
        if opts.nogaps and a.is_gap:
            continue
        if opts.bed12:
            line = a.bedline
            a.object_beg, a.object_end = ts, te
            line += "\t" + a.bedextra
            print(line, file=fw)
        elif opts.gff:
            print(a.gffline(), file=fw)
        else:
            print(a.bedline, file=fw)


def gaps(args):
    """
    %prog gaps agpfile

    Print out the distribution of gapsizes. Option --merge allows merging of
    adjacent gaps which is used by tidy().
    """
    from jcvi.graphics.histogram import loghistogram

    p = OptionParser(gaps.__doc__)
    p.add_option(
        "--merge",
        dest="merge",
        default=False,
        action="store_true",
        help="Merge adjacent gaps (to conform to AGP specification)",
    )
    p.add_option(
        "--header",
        default=False,
        action="store_true",
        help="Produce an AGP header",
    )

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    merge = opts.merge
    (agpfile,) = args

    if merge:
        merged_agpfile = agpfile.replace(".agp", ".merged.agp")
        fw = open(merged_agpfile, "w")

    agp = AGP(agpfile)
    sizes = []
    data = []  # store merged AGPLine's
    priorities = ("centromere", "telomere", "scaffold", "contig", "clone", "fragment")

    for is_gap, alines in groupby(agp, key=lambda x: (x.object, x.is_gap)):
        alines = list(alines)
        is_gap = is_gap[1]
        if is_gap:
            gap_size = sum(x.gap_length for x in alines)
            gap_types = set(x.gap_type for x in alines)
            for gtype in ("centromere", "telomere"):
                if gtype in gap_types:
                    gap_size = gtype

            sizes.append(gap_size)
            b = deepcopy(alines[0])
            b.object_beg = min(x.object_beg for x in alines)
            b.object_end = max(x.object_end for x in alines)
            b.gap_length = sum(x.gap_length for x in alines)

            assert b.gap_length == b.object_end - b.object_beg + 1
            b.component_type = "U" if b.gap_length == 100 else "N"

            gtypes = [x.gap_type for x in alines]
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

    loghistogram(sizes)

    if opts.header:
        AGP.print_header(
            fw,
            organism="Medicago truncatula",
            taxid=3880,
            source="J. Craig Venter Institute",
        )

    if merge:
        for ob, bb in groupby(data, lambda x: x.object):
            for i, b in enumerate(bb):
                b.part_number = i + 1
                print(b, file=fw)
        return merged_agpfile


def tidy(args):
    """
    %prog tidy agpfile componentfasta

    Given an agp file, run through the following steps:
    1. Trim components with dangling N's
    2. Merge adjacent gaps
    3. Trim gaps at the end of an object
    4. Reindex the agp

    Final output is in `.tidy.agp`.
    """
    p = OptionParser(tidy.__doc__)
    p.add_option(
        "--nogaps",
        default=False,
        action="store_true",
        help="Remove all gap lines",
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    agpfile, componentfasta = args
    originalagpfile = agpfile

    # Step 1: Trim terminal Ns
    tmpfasta = "tmp.fasta"
    trimmed_agpfile = build(
        [agpfile, componentfasta, tmpfasta, "--newagp", "--novalidate"]
    )
    os.remove(tmpfasta)
    agpfile = trimmed_agpfile
    agpfile = reindex([agpfile, "--inplace"])

    # Step 2: Merge adjacent gaps
    merged_agpfile = gaps([agpfile, "--merge"])
    os.remove(agpfile)

    # Step 3: Trim gaps at the end of object
    agpfile = merged_agpfile
    agp = AGP(agpfile)
    newagpfile = agpfile.replace(".agp", ".fixed.agp")
    fw = open(newagpfile, "w")
    for object, a in groupby(agp, key=lambda x: x.object):
        a = list(a)
        if a[0].is_gap:
            g, a = a[0], a[1:]
            logging.debug("Trim beginning Ns({0}) of {1}".format(g.gap_length, object))
        if a and a[-1].is_gap:
            a, g = a[:-1], a[-1]
            logging.debug("Trim trailing Ns({0}) of {1}".format(g.gap_length, object))
        print("\n".join(str(x) for x in a), file=fw)
    fw.close()
    os.remove(agpfile)

    # Step 4: Final reindex
    agpfile = newagpfile
    reindex_opts = [agpfile, "--inplace"]
    if opts.nogaps:
        reindex_opts += ["--nogaps"]
    agpfile = reindex(reindex_opts)

    tidyagpfile = originalagpfile.replace(".agp", ".tidy.agp")
    shutil.move(agpfile, tidyagpfile)

    logging.debug("File written to `%s`.", tidyagpfile)
    return tidyagpfile


def build(args):
    """
    %prog build agpfile componentfasta targetfasta

    Build targetfasta based on info from agpfile
    """
    p = OptionParser(build.__doc__)
    p.add_option(
        "--newagp",
        dest="newagp",
        default=False,
        action="store_true",
        help="Check components to trim dangling N's",
    )
    p.add_option(
        "--novalidate",
        dest="novalidate",
        default=False,
        action="store_true",
        help="Don't validate the agpfile",
    )
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    agpfile, componentfasta, targetfasta = args
    validate = not opts.novalidate

    if opts.newagp:
        assert agpfile.endswith(".agp")
        newagpfile = agpfile.replace(".agp", ".trimmed.agp")
        newagp = open(newagpfile, "w")
    else:
        newagpfile = None
        newagp = None

    agp = AGP(agpfile, validate=validate, sorted=True)
    agp.build_all(componentfasta=componentfasta, targetfasta=targetfasta, newagp=newagp)
    logging.debug("Target fasta written to `%s`.", targetfasta)

    return newagpfile


def validate(args):
    """
    %prog validate agpfile componentfasta targetfasta

    validate consistency between agpfile and targetfasta
    """
    p = OptionParser(validate.__doc__)

    opts, args = p.parse_args(args)

    try:
        agpfile, componentfasta, targetfasta = args
    except Exception as e:
        sys.exit(p.print_help())

    agp = AGP(agpfile)
    build = Fasta(targetfasta)
    bacs = Fasta(componentfasta, index=False)

    # go through this line by line
    for aline in agp:
        try:
            build_seq = build.sequence(
                dict(chr=aline.object, start=aline.object_beg, stop=aline.object_end)
            )

            if aline.is_gap:
                assert build_seq.upper() == aline.gap_length * "N", (
                    "gap mismatch: %s" % aline
                )
            else:
                bac_seq = bacs.sequence(
                    dict(
                        chr=aline.component_id,
                        start=aline.component_beg,
                        stop=aline.component_end,
                        strand=aline.orientation,
                    )
                )

                assert build_seq.upper() == bac_seq.upper(), (
                    "sequence mismatch: %s" % aline
                )

            logging.debug(
                "%s:%d-%d verified" % (aline.object, aline.object_beg, aline.object_end)
            )

        except Exception as e:
            logging.error(e)


if __name__ == "__main__":
    main()
