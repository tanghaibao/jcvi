#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Script to write and assist SOAPdenovo assembly.
"""

import os.path as op
import sys

from itertools import groupby
from optparse import OptionParser

from jcvi.assembly.base import FastqNamings, Library
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('prepare', 'prepare SOAP config files and run script'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def prepare(args):
    """
    %prog prepare *.fastq

    Scan input fastq files (see below) and write SOAP config files based
    on inputfiles.
    """
    from jcvi.utils.iter import grouper

    p = OptionParser(prepare.__doc__ + FastqNamings)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fnames = args
    for x in fnames:
        assert op.exists(x), "File `{0}` not found.".format(x)

    cfgfile = "soap.config"
    fw = open(cfgfile, "w")

    library_name = lambda x: "-".join(\
                op.basename(x).split(".")[0].split("-")[:2])
    libs = [(Library(x), sorted(fs)) for x, fs in \
                groupby(fnames, key=library_name)]

    libs.sort(key=lambda x: x[0].size)
    rank = 0
    singletons = []
    for lib, fs in libs:
        size = lib.size
        if size == 0:
            singletons = fs
            continue

        rank += 1
        block = "[LIB]\n"
        block += "avg_ins={0}\n".format(size)
        f = fs[0]
        reverse_seq = 0 if ".corr." in f else lib.reverse_seq
        block += "reverse_seq={0}\n".format(reverse_seq)
        block += "asm_flags={0}\n".format(lib.asm_flags)
        block += "rank={0}\n".format(rank)
        if singletons:
            fs += singletons
            singletons = []

        for f in fs:
            if ".1." in f:
                tag = "q1"
            elif ".2." in f:
                tag = "q2"
            else:
                tag = "q"
            block += "{0}={1}\n".format(tag, f)
        print >> sys.stderr, block
        print >> fw, block


if __name__ == '__main__':
    main()
