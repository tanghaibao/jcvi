#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Read and write EXCEL file.

http://www.simplistix.co.uk/presentations/python-excel.pdf

Library dependency: xlutils
"""

import sys

from jcvi.apps.base import MOptionParser

from xlrd import open_workbook
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('csv', 'Convert EXCEL to csv file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def csv(args):
    """
    %prog csv excelfile

    Convert EXCEL to csv file.
    """
    p = MOptionParser(csv.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    excelfile, = args
    csvfile = excelfile.rsplit(".", 1)[0] + ".csv"
    wb = open_workbook(excelfile)
    fw = open(csvfile, "w")
    for s in wb.sheets():
        print >> sys.stderr, 'Sheet:',s.name
        for row in range(s.nrows):
            values = []
            for col in range(s.ncols):
                values.append(s.cell(row, col).value)
            print >> fw, ','.join(str(x) for x in values)


if __name__ == '__main__':
    main()
