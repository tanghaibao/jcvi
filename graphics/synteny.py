#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog mcscan.txt all.bed

Illustrate MCscan multiple collinearity alignments.
"""


import sys
import logging

from optparse import OptionParser

from jcvi.algorithms.synteny import BlockFile
from jcvi.formats.bed import Bed
from jcvi.formats.base import LineFile
from jcvi.graphics.base import plt, _, set_image_options
from jcvi.apps.base import debug
debug()


class Layout (LineFile):

    def __init__(self, filename):
        super(Layout, self).__init__(filename)


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = set_image_options(p)

    if len(args) != 2:
        sys.exit(not p.print_help())

    datafile, bedfile = args
    bed = Bed(bedfile)
    order = bed.order
    bf = BlockFile(datafile)
    for i in xrange(bf.ncols):
        bf.get_extent(i, order)

    return

    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


if __name__ == '__main__':
    main()
