#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Genbank AGP file format, see spec here
http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml
"""

import sys
import itertools
import logging

from optparse import OptionParser

from jcvi.formats.base import LineFile
from jcvi.formats.fasta import Fasta
from jcvi.utils.iter import pairwise
from jcvi.apps.base import ActionDispatcher, set_debug


Valid_component_type = "ADFGNOPUW"
Valid_gap_type = ("fragment", "clone", "contig", "centromere", "short_arm",
        "heterochromatin", "telomere", "repeat")
Valid_orientation = ("+", "-", "0", "na")


class AGPLine (object):

    def __init__(self, row, validate=True):

        atoms = row.split('\t')
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
        # bed formatted line
        gid = self.component_id if not self.is_gap else "gap"
        return "\t".join((self.object, str(self.object_beg-1),
                str(self.object_end), gid, '1000', self.orientation))

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

        self.sort(key=lambda x: (x.object, x.object_beg))
        
        if validate:
            self.validate_all()

    @property
    def simple_agp(self):
        """
        returns a dict with component_id => agpline
        """
        return dict((x.component_id, x) for x in self if not x.is_gap)


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
        for ob, lines_with_same_ob in itertools.groupby(self, 
                key=lambda x: x.object):
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

            assert total_bp == line.object_end, \
                    "cumulative base pairs (%d) does not match (%d)" % \
                    (total_bp, line.object_end)


        print >>fw, ">%s\n%s" % (object, ''.join(components))
        logging.debug("Write object %s to fasta %s" % (object, fw.name))


    def build_all(self, componentfasta, targetfasta):

        from jcvi.formats.fasta import Fasta

        f = Fasta(componentfasta, index=False)
        fw = open(targetfasta, "w")

        for ob, lines_with_same_ob in itertools.groupby(self, 
                key=lambda x: x.object):

            lines = list(lines_with_same_ob)
            self.build_one(ob, lines, f, fw)


def main():

    actions = (
        ('build', 'given agp file and component fasta file, build the' + \
                 'pseudomolecule fasta'),
        ('validate', 'given agp file, component fasta and pseudomolecule fasta, ' + \
                     'validate if the build is correct')
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def build(args):
    """
    %prog build agpfile componentfasta targetfasta
    
    build targetfasta based on info from agpfile
    """
    p = OptionParser(build.__doc__)

    set_debug(p, args)
    opts, args = p.parse_args(args)

    try:
        agpfile, componentfasta, targetfasta  = args
    except Exception, e:
        logging.error(str(e))
        sys.exit(p.print_help())

    agp = AGP(agpfile)
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
        logging.error(str(e))
        sys.exit(p.print_help())

    agp = AGP(agpfile)
    build = Fasta(targetfasta, key_function=lambda x: x[0], index=False)
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
