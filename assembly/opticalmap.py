#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Optical map alignment parser.
"""

import sys

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

    def write_bed(self, fw=sys.stdout):
        for a in self.alignments:
            reference_map_name = a.reference_map_name
            aligned_map_name = a.aligned_map_name

            ref_map = self.maps[reference_map_name]
            aligned_map = self.maps[aligned_map_name]

            ref_blocks = ref_map.cumsizes
            aligned_blocks = aligned_map.cumsizes

            score = a.soma_score
            orientation = a.orientation
            for i, l, r in a.alignment:
                start = 0 if i == 0 else (aligned_blocks[i - 1] - 1)
                end = aligned_blocks[i] - 1
                accn = "{0}:{1}-{2}".format(aligned_map_name,
                        start, end)

                start = ref_blocks[l - 1] - 1
                end = ref_blocks[r] - 1
                print >> fw, "\t".join(str(x) for x in \
                        (reference_map_name, start, end,
                         accn, "{0:.1f}".format(score), orientation))


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
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed xmlfile

    Print summary of optical map alignment in BED format.
    """
    from jcvi.formats.bed import sort

    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    xmlfile, = args
    bedfile = xmlfile.rsplit(".", 1)[0] + ".bed"
    fw = open(bedfile, "w")

    om = OpticalMap(xmlfile)
    om.write_bed(fw)
    fw.close()

    sort([bedfile, "--inplace"])


if __name__ == '__main__':
    main()
