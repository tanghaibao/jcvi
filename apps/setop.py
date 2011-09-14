#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Set operations on files.
"""

import sys

from optparse import OptionParser

from jcvi.formats.base import BaseFile
from jcvi.apps.base import ActionDispatcher, debug
debug()


class SetFile (BaseFile, set):

    def __init__(self, filename):
        super(SetFile, self).__init__(filename)
        fp = open(filename)
        for x in fp:
            self.add(x.split()[0])


def main():

    actions = (
        ('setop', 'set operations on files'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    statement, = args
    fa, op, fb = statement.split()
    assert op in ('|', '&', '-', '^')

    fa = SetFile(fa)
    fb = SetFile(fb)

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
