#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# paf.py
# formats
#
# Created by Haibao Tang on 09/03/20
# Copyright © 2020 Haibao Tang. All rights reserved.
#

import sys
import logging

from jcvi.formats.base import must_open
from jcvi.apps.base import OptionParser, ActionDispatcher


class PAFLine:
    """
    PAF specification
    https://github.com/lh3/miniasm/blob/master/PAF.md
    """

    __slots__ = (
        "query",
        "qsize",
        "qstart",
        "qstop",
        "orientation",
        "subject",
        "ssize",
        "sstart",
        "sstop",
        "nmatch",
        "hitlen",
        "mapq",
    )

    def __init__(self, row):
        args = row.split()
        self.query = args[0]
        self.qsize = int(args[1])
        self.qstart = int(args[2]) + 1
        self.qstop = int(args[3])
        self.orientation = args[4]
        self.subject = args[5]
        self.ssize = int(args[6])
        self.sstart = int(args[7]) + 1
        self.sstop = int(args[8])
        self.nmatch = int(args[9])
        self.hitlen = int(args[10])
        self.mapq = int(args[11])

    @property
    def sbedline(self):
        return "\t".join(
            str(x)
            for x in (
                self.subject,
                self.sstart - 1,
                self.sstop,
                self.query,
                self.hitlen,
                self.orientation,
            )
        )

    @property
    def qbedline(self):
        return "\t".join(
            str(x)
            for x in (
                self.query,
                self.qstart - 1,
                self.qstop,
                self.subject,
                self.hitlen,
                self.orientation,
            )
        )


def bed(args):
    """
    %prog bed paffile

    Print out BED file based on coordinates in BLAST PAF results. By default,
    write out subject positions. Use --swap to write query positions.
    """
    from jcvi.formats.bed import sort as sort_bed

    p = OptionParser(bed.__doc__)
    p.add_option(
        "--swap", default=False, action="store_true", help="Write query positions"
    )

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (paffile,) = args
    write_qbed = opts.swap
    bedfile = "{}.{}.bed".format(
        paffile.rsplit(".", 1)[0], "query" if write_qbed else "subject"
    )
    with must_open(paffile) as fp, open(bedfile, "w") as fw:
        for row in fp:
            b = PAFLine(row)
            if write_qbed:
                print(b.qbedline, file=fw)
            else:
                print(b.sbedline, file=fw)

    logging.debug("File written to `%s`.", bedfile)
    sort_bed([bedfile, "-i"])
    return bedfile


def main():
    actions = (("bed", "get BED file from PAF"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
