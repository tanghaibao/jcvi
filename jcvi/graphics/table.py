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

from wand.image import Image

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
            reader = csv.reader(csvfile, skipinitialspace=True)
            self.header = [markup(x) for x in next(reader)]
            self.append(self.header)
            for row in reader:
                is_image_file = row[0].startswith("file://")
                if is_image_file:
                    images = []
                    for filenames in row:
                        images.append([Image(filename=f) for f in filenames.split("|")])
                    self.append(images)
                else:
                    self.append(row)
        print(self.header)
        print(self)

    @property
    def rows(self):
        return len(self)

    @property
    def columns(self):
        return len(self.header)


def draw_table(ax, csv_table):
    rows = csv_table.rows
    columns = csv_table.columns
    xinterval = 1.0 / columns
    yinterval = 1.0 / rows
    for i, row in enumerate(csv_table):
        for j, cell in enumerate(row):
            if isinstance(cell, list):
                continue
            else:
                ax.text(
                    (j + 0.5) * xinterval,
                    1 - (i + 0.5) * yinterval,
                    cell,
                    ha="center",
                    va="center",
                )


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
    pf = csvfile.rsplit(".", 1)[0]

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    csv_table = CsvTable(csvfile)

    draw_table(root, csv_table)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main(sys.argv[1:])
