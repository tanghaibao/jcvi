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

    def column_widths(self, total=1):
        # Get the maximum width for each column
        max_widths = [0] * self.columns
        for row in self:
            for j, cell in enumerate(row):
                if isinstance(cell, list):
                    continue
                max_widths[j] = max(max_widths[j], len(cell))
        total_width = sum(max_widths)
        return [x * total / total_width for x in max_widths]

    @property
    def rows(self):
        return len(self)

    @property
    def columns(self):
        return len(self.header)


def draw_multiple_images_in_rectangle(ax, images, rect, box_width, yinflation=1):
    """ Draw multiple images in given rectangle. Used by draw_table().

    Args:
        ax (matplotlib axes): matplotlib axes
        images (List[image]): List of images
        rect (Tuple[float]): (left, bottom, width, height)
        box_width (float): Width of the image square
    """
    n_images = len(images)
    left, bottom, width, height = rect
    box_start = (width - n_images * box_width) / 2
    left += box_start
    bottom += (height - box_width * yinflation) / 2
    for image in images:
        extent = (left, left + box_width, bottom, bottom + box_width * yinflation)
        ax.imshow(image, extent=extent, aspect="auto")
        left += box_width


def draw_table(ax, csv_table, extent=(0, 1, 0, 1), stripe_color="beige", yinflation=1):
    """ Draw table on canvas.

    Args:
        ax (matplotlib axes): matplotlib axes
        csv_table (CsvTable): Parsed CSV table
        extent (tuple, optional): (left, right, bottom, top). Defaults to (0, 1, 0, 1).
        stripe_color (str, optional): Stripe color of the table. Defaults to
        "beige".
        yinflation (float, optional): Inflate on y since imshow aspect ratio
        sometimes create warped images. Defaults to 1.
    """
    left, right, bottom, top = extent
    width = right - left
    height = top - bottom
    rows = csv_table.rows
    columns = csv_table.columns
    column_widths = csv_table.column_widths(width)
    print(column_widths)

    xinterval = width / columns
    yinterval = height / rows
    for i, row in enumerate(csv_table):
        should_stripe = i % 2 == 0
        contain_images = isinstance(row[0], list)
        xstart = left
        if contain_images:
            box_width = min(
                min(column_widths[j] / len(x) for j, x in enumerate(row)), yinterval
            )
        for j, cell in enumerate(row):
            xinterval = column_widths[j]
            xmid = xstart + xinterval / 2
            ymid = top - (i + 0.5) * yinterval
            if contain_images:
                # There may be multiple images, center them
                rect = (xstart, top - (i + 1) * yinterval, xinterval, yinterval)
                draw_multiple_images_in_rectangle(
                    ax, cell, rect, box_width, yinflation=yinflation
                )
                should_stripe = False
            else:
                ax.text(
                    xmid, ymid, cell, ha="center", va="center",
                )

            xstart += column_widths[j]

        if not should_stripe:
            continue

        # Draw the stripes, extend a little longer horizontally
        xpad = 0.01
        ax.add_patch(
            Rectangle(
                (left - xpad, top - (i + 1) * yinterval),
                width + 2 * xpad,
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

    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main(sys.argv[1:])
