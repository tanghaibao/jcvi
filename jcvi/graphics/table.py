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
    load_image,
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
                        images.append(
                            [
                                load_image(filename=f.replace("file://", ""))
                                for f in filenames.split("|")
                            ]
                        )
                    self.append(images)
                else:
                    self.append(row)
        print(self.header)
        # print(self)

    @property
    def rows(self):
        return len(self)

    @property
    def columns(self):
        return len(self.header)


def draw_multiple_images_in_rectangle(ax, images, rect):
    n_images = len(images)
    left, bottom, width, height = rect
    box_width = min(width / n_images, height)
    box_start = (width - n_images * box_width) / 2
    left += box_start
    bottom += (height - box_width) / 2
    for image in images:
        extent = (left, left + box_width, bottom, bottom + box_width)
        ax.imshow(image, extent=extent)
        left += box_width
        print(extent)


def draw_table(ax, csv_table, stripe_color="beige"):
    rows = csv_table.rows
    columns = csv_table.columns
    xinterval = 1.0 / columns
    yinterval = 1.0 / rows
    for i, row in enumerate(csv_table):
        should_stripe = i % 2 == 0
        for j, cell in enumerate(row):
            xmid = (j + 0.5) * xinterval
            ymid = 1 - (i + 0.5) * yinterval
            if isinstance(cell, list):
                # There may be multiple images, center them
                rect = (j * xinterval, 1 - (i + 1) * yinterval, xinterval, yinterval)
                print(rect)
                draw_multiple_images_in_rectangle(ax, cell, rect)
                should_stripe = False
            else:
                ax.text(
                    xmid, ymid, cell, ha="center", va="center",
                )

        if not should_stripe:
            continue

        # Draw the stripes
        ax.add_patch(
            Rectangle(
                (0, 1 - (i + 1) * yinterval),
                1,
                yinterval,
                fc=stripe_color,
                ec=stripe_color,
            )
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
