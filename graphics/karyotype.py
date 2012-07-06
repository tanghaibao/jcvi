#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog anchorsfile tracksfile all.bed

Illustrate macrosynteny between tracks which represent individual genomes.
"""


import sys
import logging

from optparse import OptionParser

from jcvi.graphics.base import plt, _, set_image_options
from jcvi.apps.base import debug
debug()


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = set_image_options(p)

    if len(args) != 3:
        sys.exit(not p.print_help())

    datafile, = args
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
