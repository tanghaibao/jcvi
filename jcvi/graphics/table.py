#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# table.py
# graphics
#
# Created by Haibao Tang on 05/25/20
# Copyright © 2020 Haibao Tang. All rights reserved.
#
import csv
import sys
import os.path as op

from jcvi.apps.base import OptionParser
from jcvi.graphics.base import (
    FancyBboxPatch,
    Rectangle,
    linear_shade,
    markup,
    normalize_axes,
    plt,
    savefig,
)


class CsvTable(list):
    def __init__(self, csvfile="table.csv"):
        with open(csvfile) as csvfile:
            reader = csv.reader(csvfile)


def main(args):
    """
    %prog table.csv

    Render a table on canvas. Input is a CSV file.
    """
    p = OptionParser(main.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="7x7")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (csvfile,) = args

    pf = op.basename(csvfile).split(".")[0]
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main(sys.argv[1:])
