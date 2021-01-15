#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Read and write EXCEL file.

http://www.simplistix.co.uk/presentations/python-excel.pdf

Library dependency: xlutils
"""
import os.path as op
import sys
import logging

from jcvi.apps.base import OptionParser, ActionDispatcher


class ColorMatcher(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self.unused_colors = set(self.xlwt_colors)
        # Never use black.
        self.unused_colors.discard((0, 0, 0))

    # Culled from a table at http://www.mvps.org/dmcritchie/excel/colors.htm
    xlwt_colors = [
        (0, 0, 0),
        (255, 255, 255),
        (255, 0, 0),
        (0, 255, 0),
        (0, 0, 255),
        (255, 255, 0),
        (255, 0, 255),
        (0, 255, 255),
        (0, 0, 0),
        (255, 255, 255),
        (255, 0, 0),
        (0, 255, 0),
        (0, 0, 255),
        (255, 255, 0),
        (255, 0, 255),
        (0, 255, 255),
        (128, 0, 0),
        (0, 128, 0),
        (0, 0, 128),
        (128, 128, 0),
        (128, 0, 128),
        (0, 128, 128),
        (192, 192, 192),
        (128, 128, 128),
        (153, 153, 255),
        (153, 51, 102),
        (255, 255, 204),
        (204, 255, 255),
        (102, 0, 102),
        (255, 128, 128),
        (0, 102, 204),
        (204, 204, 255),
        (0, 0, 128),
        (255, 0, 255),
        (255, 255, 0),
        (0, 255, 255),
        (128, 0, 128),
        (128, 0, 0),
        (0, 128, 128),
        (0, 0, 255),
        (0, 204, 255),
        (204, 255, 255),
        (204, 255, 204),
        (255, 255, 153),
        (153, 204, 255),
        (255, 153, 204),
        (204, 153, 255),
        (255, 204, 153),
        (51, 102, 255),
        (51, 204, 204),
        (153, 204, 0),
        (255, 204, 0),
        (255, 153, 0),
        (255, 102, 0),
        (102, 102, 153),
        (150, 150, 150),
        (0, 51, 102),
        (51, 153, 102),
        (0, 51, 0),
        (51, 51, 0),
        (153, 51, 0),
        (153, 51, 102),
        (51, 51, 153),
        (51, 51, 51),
    ]

    @staticmethod
    def color_distance(rgb1, rgb2):
        # Adapted from Colour metric by Thiadmer Riemersma,
        # http://www.compuphase.com/cmetric.htm
        rmean = (rgb1[0] + rgb2[0]) / 2
        r = rgb1[0] - rgb2[0]
        g = rgb1[1] - rgb2[1]
        b = rgb1[2] - rgb2[2]
        return (
            (((512 + rmean) * r * r) / 256)
            + 4 * g * g
            + (((767 - rmean) * b * b) / 256)
        )

    def match_color_index(self, color):
        """Takes an "R,G,B" string or wx.Color and returns a matching xlwt
        color.
        """
        from jcvi.utils.webcolors import color_diff

        if isinstance(color, int):
            return color
        if color:
            if isinstance(color, str):
                rgb = map(int, color.split(","))
            else:
                rgb = color.Get()
            logging.disable(logging.DEBUG)
            distances = [color_diff(rgb, x) for x in self.xlwt_colors]
            logging.disable(logging.NOTSET)
            result = distances.index(min(distances))
            self.unused_colors.discard(self.xlwt_colors[result])
            return result

    def get_unused_color(self):
        """Returns an xlwt color index that has not been previously returned by
        this instance.  Attempts to maximize the distance between the color and
        all previously used colors.
        """
        if not self.unused_colors:
            # If we somehow run out of colors, reset the color matcher.
            self.reset()
        used_colors = [c for c in self.xlwt_colors if c not in self.unused_colors]
        result_color = max(
            self.unused_colors,
            key=lambda c: min(self.color_distance(c, c2) for c2 in used_colors),
        )
        result_index = self.xlwt_colors.index(result_color)
        self.unused_colors.discard(result_color)
        return result_index


def main():

    actions = (
        ("csv", "Convert EXCEL to csv file"),
        ("fromcsv", "Convert csv file to EXCEL"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fromcsv(args):
    """
    %prog fromcsv csvfile

    Convert csv file to EXCEL.
    """
    from csv import reader
    from xlwt import Workbook, easyxf
    from jcvi.formats.base import flexible_cast

    p = OptionParser(fromcsv.__doc__)
    p.add_option(
        "--noheader",
        default=False,
        action="store_true",
        help="Do not treat the first row as header",
    )
    p.add_option("--rgb", default=-1, type="int", help="Show RGB color box")
    p.set_sep()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (csvfile,) = args
    header = not opts.noheader
    rgb = opts.rgb
    excelfile = csvfile.rsplit(".", 1)[0] + ".xls"

    data = []
    for row in reader(open(csvfile), delimiter=opts.sep):
        data.append(row)

    w = Workbook()
    s = w.add_sheet(op.basename(csvfile))

    header_style = easyxf("font: bold on")
    if header:
        s.panes_frozen = True
        s.horz_split_pos = 1

    cm = ColorMatcher()
    for i, row in enumerate(data):
        for j, cell in enumerate(row):
            cell = flexible_cast(cell)
            if header and i == 0:
                s.write(i, j, cell, header_style)
            else:
                if j == rgb:
                    cix = cm.match_color_index(cell)
                    color_style = easyxf("font: color_index {0}".format(cix))
                    s.write(i, j, cell, color_style)
                else:
                    s.write(i, j, cell)

    w.save(excelfile)
    logging.debug("File written to `%s`.", excelfile)
    return excelfile


def csv(args):
    """
    %prog csv excelfile

    Convert EXCEL to csv file.
    """
    from xlrd import open_workbook

    p = OptionParser(csv.__doc__)
    p.set_sep(sep=",")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (excelfile,) = args
    sep = opts.sep
    csvfile = excelfile.rsplit(".", 1)[0] + ".csv"
    wb = open_workbook(excelfile)
    fw = open(csvfile, "w")
    for s in wb.sheets():
        print("Sheet:", s.name, file=sys.stderr)
        for row in range(s.nrows):
            values = []
            for col in range(s.ncols):
                values.append(s.cell(row, col).value)
            print(sep.join(str(x) for x in values), file=fw)


if __name__ == "__main__":
    main()
