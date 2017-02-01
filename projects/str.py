#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Related scripts for the HLI-STR (TREDPARSE) paper.
"""


import sys
import logging

from jcvi.graphics.base import plt, savefig, normalize_axes
from jcvi.apps.base import OptionParser


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = p.set_image_options()

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
