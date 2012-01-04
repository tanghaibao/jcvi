#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Set operations on files.
"""

import sys

from optparse import OptionParser

from jcvi.formats.base import BaseFile, DictFile, must_open
from jcvi.apps.base import ActionDispatcher, debug, set_outfile
debug()


class SetFile (BaseFile, set):

    def __init__(self, filename, column=0):
        super(SetFile, self).__init__(filename)
        fp = open(filename)
        for x in fp:
            key = x.split()[column] if column >= 0 else x.strip()
            self.add(key)


def main():

    actions = (
        ('setop', 'set operations on files'),
        ('join', 'join tabular files based on common column'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def join(args):
    """
    %prog join file1.txt file2.txt ..

    Join tabular files based on common column.
    """
    from jcvi.utils.iter import flatten

    p = OptionParser(join.__doc__)
    p.add_option("--column", default=0, type="int",
                 help="The column to pivot on, 0-based [default: %default]")
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    c = opts.column

    # Maintain the first file line order, and combine other files into it
    pivotfile = args[0]
    files = [DictFile(f, keypos=c, valuepos=None, delimiter="\t") \
                        for f in args]
    otherfiles = files[1:]
    header = "\t".join(flatten([x.filename] * x.ncols for x in files))

    fp = open(pivotfile)
    fw = must_open(opts.outfile, "w")
    print >> fw, header
    for row in fp:
        row = row.rstrip()
        atoms = row.split("\t")
        newrow = atoms
        key = atoms[c]
        for d in otherfiles:
            drow = d.get(key, ["na"] * d.ncols)
            newrow += drow
        print >> fw, "\t".join(newrow)


def setop(args):
    """
    %prog setop "fileA & fileB" > newfile

    Perform set operations, except on files. The files (fileA and fileB) contain
    list of ids. The operator is one of the four:

    |: union (elements found in either file)
    &: intersection (elements found in both)
    -: difference (elements in fileA but not in fileB)
    ^: symmetric difference (elementes found in either set but not both)

    Please quote the argument to avoid shell interpreting | and &.
    """
    p = OptionParser(setop.__doc__)
    p.add_option("--column", default=0, type="int",
                 help="The column to extract, 0-based, -1 to disable [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    statement, = args
    fa, op, fb = statement.split()
    assert op in ('|', '&', '-', '^')

    column = opts.column
    fa = SetFile(fa, column=column)
    fb = SetFile(fb, column=column)

    if op == '|':
        t = fa | fb
    elif op == '&':
        t = fa & fb
    elif op == '-':
        t = fa - fb
    elif op == '^':
        t = fa ^ fb

    for x in sorted(t):
        print x


if __name__ == '__main__':
    main()
