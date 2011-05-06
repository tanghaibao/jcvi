#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog blastfile fastafiles

Calculate the vector clear range file based on the blastfile to the vectors
"""

import sys

from optparse import OptionParser

from jcvi.formats.blast import Blast
from jcvi.formats.fasta import Fasta
from jcvi.utils.range import range_minmax
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 0:
        sys.exit(p.print_help())

    blastfile = args[0]
    fastafiles = args[1:]

    sizes = {}
    for fa in fastafiles:
        f = Fasta(fa)
        sizes.update(f.itersizes())

    b = Blast(blastfile)
    seen = set()
    for query, hits in b.iter_hits():
        #if query not in sizes: continue

        qsize = sizes[query]
        vectors = list((x.qstart, x.qstop) for x in hits)
        vmin, vmax = range_minmax(vectors)

        left_size = vmin - 1
        right_size = qsize - vmax

        if left_size > right_size:
            clr_start, clr_end = 0, vmin
        else:
            clr_start, clr_end = vmax, qsize

        print "\t".join(str(x) for x in (query, clr_start, clr_end))
        del sizes[query]

    for q, size in sorted(sizes.items()):
        print "\t".join(str(x) for x in (q, 0, size))


if __name__ == '__main__':
    main()
