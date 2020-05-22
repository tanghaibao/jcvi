#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# phylogeny.py
# graphics
#
# Created by Haibao Tang on 05/21/20
# Copyright © 2020 Haibao Tang. All rights reserved.
#
import csv
import sys

from jcvi.apps.base import ActionDispatcher, OptionParser


def lcn(args):
    """
    %prog lcn Orthogroups/Orthogroups.tsv
    """
    p = OptionParser(lcn.__doc__)
    p.add_option("--min-single-ratio", default=0.9, help="Single copy ratio must be > ")
    p.add_option("--max-zero-ratio", default=0.1, help="Zero copy ratio must be < ")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (groups_tsv,) = args
    total = 0
    with open(groups_tsv) as fp:
        reader = csv.reader(fp, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            counts = [len(x.split(", ")) if x.strip() != "" else 0 for x in row[1:]]
            single_ratio = sum([x == 1 for x in counts]) / len(counts)
            zero_ratio = sum([x == 0 for x in counts]) / len(counts)
            if single_ratio <= opts.min_single_ratio:
                continue
            if zero_ratio >= opts.max_zero_ratio:
                continue
            print(row[0], single_ratio, zero_ratio, counts)
            total += 1
    print(total)


def main():
    actions = (("lcn", "collect low copy ortholog groups from OrthoFinder results"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
