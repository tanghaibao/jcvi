#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
From genomeCovergeBed results, initialize the count array, set cutoffs
and optimize against the truth, to determine the cutoff for incorporating
RNA-seq into annotation pipelines.
"""

import sys
import os.path as op
import numpy as np
import logging

from itertools import groupby
from optparse import OptionParser

from jcvi.formats.sizes import Sizes
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('count', 'initialize the count array'),
        ('query', 'query the count array to get depth at particular site'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def query(args):
    """
    %prog query binfile fastafile ctgID baseID

    Get the depth at a particular base.
    """
    p = OptionParser(query.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    countsfile, fastafile, ctgID, baseID = args
    if not op.exists(countsfile):
        logging.error("Binary file `{0}` not found. Rerun count()."\
                        .format(countsfile))
        return

    fastasize, sizes, offsets = get_offsets(fastafile)
    ar = np.memmap(countsfile, dtype=np.uint8, mode="r")
    oi = offsets[ctgID] + int(baseID) - 1
    print "\t".join((ctgID, baseID, str(ar[oi])))


def update_array(ar, coveragefile, sizes, offsets):
    fp = open(coveragefile)
    logging.debug("Parse file `{0}`".format(coveragefile))
    for k, rows in groupby(fp, key=(lambda x: x.split()[0])):
        rows = list(rows)
        offset = offsets[k]
        print k, offset
        ctglen = len(rows)

        if ctglen < 100000:
            sys.stdout.write(".")

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
        logging.error("`{0}` file exists. Remove before proceed.")
        return

    fastasize, sizes, offsets = get_offsets(fastafile)
    logging.debug("Initialize array of uint8 with size {0}".format(fastasize))
    ar = np.zeros(fastasize, dtype=np.uint8)

    update_array(ar, coveragefile, sizes, offsets)

    ar.tofile(countsfile)
    logging.debug("Array written to `{0}`".format(countsfile))


if __name__ == '__main__':
    main()
