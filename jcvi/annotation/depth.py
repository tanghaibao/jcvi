#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
From genomeCovergeBed results, initialize the count array, set cutoffs
and optimize against the truth, to determine the cutoff for incorporating
RNA-seq into annotation pipelines.
"""
from __future__ import print_function

import sys
import os.path as op
import numpy as np
import logging

from itertools import groupby

from jcvi.formats.sizes import Sizes
from jcvi.formats.base import BaseFile, must_open
from jcvi.apps.base import OptionParser, ActionDispatcher


class BinFile (BaseFile):
    """
    The binfile contains per base count, fastafile provides the coordinate
    system.
    """
    def __init__(self, binfile, fastafile=None, dtype=np.uint8):
        super(BinFile, self).__init__(binfile)
        assert op.exists(binfile), \
            "Binary file `{0}` not found. Rerun depth.count().".format(binfile)
        self.dtype = dtype

    @property
    def array(self):
        binfile = self.filename
        return np.fromfile(binfile, dtype=self.dtype)

    @property
    def mmarray(self):
        binfile = self.filename
        return np.memmap(binfile, dtype=self.dtype, mode="r")


def main():

    actions = (
        ('count', 'initialize the count array'),
        ('query', 'query the count array to get depth at particular site'),
        ('merge', 'merge several count arrays into one'),
        ('bed', 'write bed files where the bases have at least certain depth',)
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed binfile fastafile

    Write bed files where the bases have at least certain depth.
    """
    p = OptionParser(bed.__doc__)
    p.add_option("-o", dest="output", default="stdout",
            help="Output file name [default: %default]")
    p.add_option("--cutoff", dest="cutoff", default=10, type="int",
            help="Minimum read depth to report intervals [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    binfile, fastafile = args
    fw = must_open(opts.output, "w")
    cutoff = opts.cutoff
    assert cutoff >= 0, "Need non-negative cutoff"

    b = BinFile(binfile)
    ar = b.array

    fastasize, sizes, offsets = get_offsets(fastafile)
    s = Sizes(fastafile)
    for ctg, ctglen in s.iter_sizes():
        offset = offsets[ctg]
        subarray = ar[offset:offset + ctglen]
        key = lambda x: x[1] >= cutoff
        for tf, array_elements in groupby(enumerate(subarray), key=key):
            array_elements = list(array_elements)
            if not tf:
                continue

            # 0-based system => 1-based system
            start = array_elements[0][0] + 1
            end = array_elements[-1][0] + 1

            mean_depth = sum([x[1] for x in array_elements]) / \
                len(array_elements)
            mean_depth = int(mean_depth)

            name = "na"
            print("\t".join(str(x) for x in (ctg, \
                    start - 1, end, name, mean_depth)), file=fw)


def merge(args):
    """
    %prog merge *.bin merged.bin

    Merge several count arrays into one. Overflows will be capped at uint8_max
    (255).
    """
    p = OptionParser(merge.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    binfiles = args[:-1]
    mergedbin = args[-1]
    if op.exists(mergedbin):
        logging.error("`{0}` file exists. Remove before proceed."\
                .format(mergedbin))
        return

    b = BinFile(binfiles[0])
    ar = b.mmarray
    fastasize, = ar.shape
    logging.debug("Initialize array of uint16 with size {0}".format(fastasize))

    merged_ar = np.zeros(fastasize, dtype=np.uint16)
    for binfile in binfiles:
        b = BinFile(binfile)
        merged_ar += b.array

    logging.debug("Resetting the count max to 255.")
    merged_ar[merged_ar > 255] = 255

    logging.debug("Compact array back to uint8 with size {0}".format(fastasize))
    merged_ar = np.array(merged_ar, dtype=np.uint8)
    merged_ar.tofile(mergedbin)
    logging.debug("Merged array written to `{0}`".format(mergedbin))


def query(args):
    """
    %prog query binfile fastafile ctgID baseID

    Get the depth at a particular base.
    """
    p = OptionParser(query.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    binfile, fastafile, ctgID, baseID = args
    b = BinFile(binfile, fastafile)
    ar = b.mmarray

    fastasize, sizes, offsets = get_offsets(fastafile)
    oi = offsets[ctgID] + int(baseID) - 1
    print("\t".join((ctgID, baseID, str(ar[oi]))))


def update_array(ar, coveragefile, sizes, offsets):
    fp = open(coveragefile)
    logging.debug("Parse file `{0}`".format(coveragefile))
    for k, rows in groupby(fp, key=(lambda x: x.split()[0])):
        rows = list(rows)
        offset = offsets[k]
        ctglen = len(rows)

        if ctglen < 100000:
            sys.stdout.write(".")
        else:
            print(k, offset)

        #assert ctglen == sizes[k]
        for i, row in enumerate(rows):
            ctgID, baseID, count = row.split()
            oi = offset + i
            newcount = ar[oi] + int(count)
            if newcount > 255:
                newcount = 255

            ar[oi] = newcount


def get_offsets(fastafile):
    s = Sizes(fastafile)
    fastasize = s.totalsize
    sizes = s.mapping
    offsets = s.cumsizes_mapping
    return fastasize, sizes, offsets


def count(args):
    """
    %prog count t.coveragePerBase fastafile

    Serialize the genomeCoverage results. The coordinate system of the count array
    will be based on the fastafile.
    """
    p = OptionParser(count.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    coveragefile, fastafile = args

    countsfile = coveragefile.split(".")[0] + ".bin"
    if op.exists(countsfile):
        logging.error("`{0}` file exists. Remove before proceed."\
                .format(countsfile))
        return

    fastasize, sizes, offsets = get_offsets(fastafile)
    logging.debug("Initialize array of uint8 with size {0}".format(fastasize))
    ar = np.zeros(fastasize, dtype=np.uint8)

    update_array(ar, coveragefile, sizes, offsets)

    ar.tofile(countsfile)
    logging.debug("Array written to `{0}`".format(countsfile))


if __name__ == '__main__':
    main()
