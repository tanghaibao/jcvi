#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
CLC bio assembly file CAS
"""

import sys

from optparse import OptionParser

from jcvi.apps.grid import GridProcess
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('split', 'split the CAS file into smaller CAS using sub_assembly'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def split(args):
    """
    %prog split casfile 1 10

    """
    p = OptionParser(split.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(p.print_help())

    casfile, start, end = args
    start = int(start)
    end = int(end)

    split_cmd = "sub_assembly -a {casfile} -o sa.{i}.cas -s {i} " + \
        "-e sa.{i}.pairs.fasta -f sa.{i}.fragments.fasta -g sa.{i}.ref.fasta"

    for i in range(start, end+1):
        cmd = split_cmd.format(casfile=casfile, i=i)
        p = GridProcess(cmd)
        p.start(path=None) # current dir


if __name__ == '__main__':
    main()
