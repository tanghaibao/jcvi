#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog fpcfile

Parser for FPC physical map file, writes a tabular file to stdout
"""
from __future__ import print_function

import sys
import logging

from jcvi.apps.base import OptionParser
from jcvi.formats.base import read_until, read_block


bac_tag = "BAC :"


class FpcRecord (object):

    def __init__(self):
        self.bac_name = ''
        self.ctg_name = ''
        self.map_left = 0
        self.map_right = 0
        self.cosmid = ''
        self.cosmid_match = ''
        self.gel_number = ''
        self.fp_number = ''
        self.bands = ''
        self.probes = []
        self.cre_date = ''
        self.mod_date = ''
        self.remark = []

    def __str__(self):
        return '\t'.join((self.bac_name, self.ctg_name,
            str(self.map_left), str(self.map_right),
            self.bands, ','.join(self.probes), ','.join(self.remark)))


class FpcReader (object):
    """
    parse fpc file, usage:

    >>> reader = FpcReader()
    >>> for rec in reader:
            print rec
    """
    def __init__(self, fpcfile):
        self._handle = open(fpcfile)

    def __iter__(self):
        line = self._handle.readline()
        if not line.startswith(bac_tag):
            read_until(self._handle, bac_tag)

        for header, seq in read_block(self._handle, bac_tag):

            rec = FpcRecord()
            assert header.startswith(bac_tag)
            rec.bac_name = header.split('\"')[1]

            for line in seq:
                if line.startswith("Map"):
                    rec.ctg_name = line.split('\"')[1]
                    if "Left" in line:
                        rec.map_left = line.split("Left")[1].split()[0]
                        rec.map_left = int(float(rec.map_left))
                    if "Right" in line:
                        rec.map_right = line.split("Right")[1].split()[0]
                        rec.map_right = int(float(rec.map_right))
                if line.startswith("Gel_number"):
                    rec.gel_number = line.split()[-1]
                if line.startswith("Fp_number"):
                    rec.fp_number = line.split()[-1]
                if line.startswith("Bands"):
                    rec.bands = line.split()[-1]
                if line.startswith("match_to_cosmid"):
                    rec.cosmid = line.split('\"')[1]
                    rec.cosmid_match = line.split("_match_")[0]
                if line.startswith("Positive"):
                    rec.probes.append(line.split('\"')[1])
                if line.startswith("Fpc_remark"):
                    rec.remark.append(line.split('\"')[1])
                if line.startswith("Creation_date"):
                    rec.cre_date = line.split('_date')[1].strip()
                if line.startswith("Modified_date"):
                    rec.mod_date = line.split('_date')[1].strip()

            yield rec


def main(fpcfile):

    fw = sys.stdout
    f = FpcReader(fpcfile)

    # first several lines are comments
    header = '\t'.join(('bac_name', 'ctg_name', 'map_left', 'map_right',
        'bands', 'probes', 'remark'))
    print(header, file=fw)

    recs = list(f)
    logging.debug("%d records parsed" % len(recs))

    recs.sort(key=lambda x: (x.ctg_name, x.map_left, x.map_right))
    for rec in recs:
        print(rec, file=fw)


if __name__ == '__main__':

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    fpcfile, = args

    main(fpcfile)
