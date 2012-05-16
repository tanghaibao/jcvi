#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Optical map alignment parser.
"""

import sys
import logging

import numpy as np

from optparse import OptionParser
from xml.etree.ElementTree import ElementTree

from jcvi.apps.base import ActionDispatcher, debug
debug()


class OpticalMap (object):

    def __init__(self, xmlfile):
        tree = ElementTree()
        self.root = tree.parse(xmlfile)
        self.maps = dict(self.iter_maps())
        self.alignments = []

        for ref, aligned, e in self.iter_alignments():
            ref_map = self.maps[ref]
            aligned_map = self.maps[aligned]
            nfrags = aligned_map.num_frags
            if e.orientation == '-':
                e.alignment = [(nfrags - i - 1, l, r) for (i, l, r) in e.alignment]
            self.alignments.append(e)

    def iter_maps(self):
        for e in self.root.findall("restriction_map"):
            e = RestrictionMap(e)
            yield e.name, e

    def iter_alignments(self):
        for e in self.root.findall("map_alignment"):
            e = MapAlignment(e)
            yield e.reference_map_name, e.aligned_map_name, e

    def write_bed(self, fw=sys.stdout, blockonly=False):
        for a in self.alignments:
            reference_map_name = a.reference_map_name
            aligned_map_name = a.aligned_map_name

            ref_map = self.maps[reference_map_name]
            aligned_map = self.maps[aligned_map_name]

            ref_blocks = ref_map.cumsizes
            aligned_blocks = aligned_map.cumsizes

            score = a.soma_score
            score = "{0:.1f}".format(score)
            orientation = a.orientation

            endpoints = []
            ref_endpoints = []
            for i, l, r in a.alignment:
                start = 0 if i == 0 else (aligned_blocks[i - 1] - 1)
                end = aligned_blocks[i] - 1
                endpoints.extend([start, end])
                accn = "{0}:{1}-{2}".format(aligned_map_name,
                        start, end)

                start = ref_blocks[l - 1] - 1
                end = ref_blocks[r] - 1
                ref_endpoints.extend([start, end])

                if not blockonly:
                    print >> fw, "\t".join(str(x) for x in \
                            (reference_map_name, start, end,
                             accn, score, orientation))

            if blockonly:
                start, end = min(endpoints), max(endpoints)
                accn = "{0}:{1}-{2}".format(aligned_map_name, start, end)

                start, end = min(ref_endpoints), max(ref_endpoints)
                print >> fw, "\t".join(str(x) for x in \
                        (reference_map_name, start, end,
                         accn, score, orientation))


class RestrictionMap (object):

    def __init__(self, node):
        num_frags = node.find("num_frags").text
        map_blocks = node.find("map_block").text

        num_frags = int(num_frags)

        self.name = node.find("name").text
        self.num_frags = num_frags
        self.map_blocks = [int(round(float(x) * 1000)) \
                           for x in map_blocks.split()]

        assert len(self.map_blocks) == self.num_frags

    @property
    def cumsizes(self):
        return np.cumsum(self.map_blocks)



class MapAlignment (object):

    def __init__(self, node):
        reference_map = node.find("reference_map")
        reference_map_name = reference_map.find("name").text

        aligned_map = node.find("aligned_map")
        aligned_map_name = aligned_map.find("name").text
        aligned_map_orientation = aligned_map.find("orientation").text

        assert aligned_map_orientation in ('N', 'R')
        self.orientation = '-' if aligned_map_orientation == 'R' else '+'

        soma_score = node.find("soma_score").text
        count = node.find("count").text

        soma_score = float(soma_score)
        count = int(count)

        self.reference_map_name = reference_map_name
        self.aligned_map_name = aligned_map_name
        self.aligned_map_orientation = aligned_map_orientation

        self.soma_score = soma_score
        self.alignment = []

        for f in node.findall("f"):
            i = f.find("i").text
            l = f.find("l").text
            r = f.find("r").text
            i, l, r = [int(x) for x in (i, l, r)]
            self.alignment.append((i, l, r))


def main():

    actions = (
        ('bed', 'convert xml format into bed format'),
        ('fasta', 'use the OM bed to scaffold and create pseudomolecules'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fasta(args):
    """
    %prog fasta bedfile scf.fasta pseudomolecules.fasta

    Use OM bed to scaffold and create pseudomolecules. bedfile can be generated
    by running jcvi.assembly.opticalmap bed --blockonly
    """
    from jcvi.formats.bed import Bed
    from jcvi.formats.sizes import Sizes
    from jcvi.formats.agp import OO, build
    from jcvi.utils.range import range_chain, range_parse, Range

    p = OptionParser(fasta.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, scffasta, pmolfasta = args
    pf = bedfile.rsplit(".", 1)[0]
    bed = Bed(bedfile)
    ranges = [Range(x.seqid, x.start, x.end, float(x.score), i) for i, x in enumerate(bed)]
    selected, score = range_chain(ranges)
    selected = [bed[x.id] for x in selected]
    oo = OO()
    seen = set()
    sizes = Sizes(scffasta).mapping
    agpfile = pf + ".agp"
    agp = open(agpfile, "w")
    for b in selected:
        scf = range_parse(b.accn).seqid
        chr = b.seqid
        cs = (chr, scf)
        if cs not in seen:
            oo.add(chr, scf, sizes[scf], b.strand)
            seen.add(cs)
        else:
            logging.debug("Seen {0}, ignored.".format(cs))

    oo.write_AGP(agp, gaptype="contig")
    build([agpfile, scffasta, pmolfasta])


def bed(args):
    """
    %prog bed xmlfile

    Print summary of optical map alignment in BED format.
    """
    from jcvi.formats.bed import sort

    p = OptionParser(bed.__doc__)
    p.add_option("--blockonly", default=False, action="store_true",
                 help="Only print out large blocks, not fragments [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    xmlfile, = args
    bedfile = xmlfile.rsplit(".", 1)[0] + ".bed"
    fw = open(bedfile, "w")

    om = OpticalMap(xmlfile)
    om.write_bed(fw, blockonly=opts.blockonly)
    fw.close()

    sort([bedfile, "--inplace"])


if __name__ == '__main__':
    main()
