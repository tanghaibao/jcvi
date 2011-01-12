#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog infile [options]

This script combines a series of pairwise alignments, sort and filter the alignments
Infile expect BLAST tabular format (-m8). 

In order to handle dups, we have to run two monotonous chains in both genomes,
first chain using ref, and a second chain using query and we will have options
to keep either the union or the intersection of retain chained alignments from
both genomes, similar to the SUPERMAP algorithm. This operation is symmetrical.
"""

import os
import sys
import collections
import logging

from optparse import OptionParser

from jcvi.formats.blast import BlastLine
from jcvi.formats.coords import CoordsLine
from jcvi.utils.range import Range, range_chain
from jcvi.apps.base import debug
debug()


def BlastOrCoordsLine(filename, filter="ref", dialect="blast"):
    assert filter in ("ref", "query")
    assert dialect in ("blast", "coords")

    filter = filter[0]
    dialect = dialect[0]
    
    fp = open(filename)
    for i, row in enumerate(fp):
        if dialect=="b":
            b = BlastLine(row)
            if filter=="q":
                query, start, end = b.query, b.qstart, b.qstop
            else:
                query, start, end = b.subject, b.sstart, b.sstop
        else:
            try:
                b = CoordsLine(row)
            except AssertionError:
                continue

            if filter=="q":
                query, start, end = b.query, b.start2, b.end2
            else:
                query, start, end = b.ref, b.start1, b.end1

        if start > end: start, end = end, start

        yield Range(query, start, end, b.score, i)


def main(blast_file, filter="intersection", dialect="blast"):
    # filter by query
    if opts.filter!="ref":
        logging.debug("filter by query")
        ranges = list(BlastOrCoordsLine(blast_file, filter="query",
            dialect=dialect))

        query_selected, query_score = range_chain(ranges)
        query_idx = set(x.id for x in query_selected)

    # filter by ref
    if filter!="query":
        logging.debug("filter by ref")
        ranges = list(BlastOrCoordsLine(blast_file, filter="ref",
            dialect=dialect))

        ref_selected, ref_score = range_chain(ranges)
        ref_idx = set(x.id for x in ref_selected)

    if filter=="ref":
        selected_idx = ref_idx

    elif filter=="query":
        selected_idx = query_idx

    elif filter=="intersection":
        logging.debug("perform intersection")
        selected_idx = ref_idx & query_idx

    elif filter=="union":
        logging.debug("perform union")
        selected_idx = ref_idx | query_idx

    assert len(selected_idx)!=0

    # selected_idx is in fact the lineno in the BLAST file
    fp = open(blast_file)
    selected_idx = iter(sorted(selected_idx))
    selected = selected_idx.next()
    for i, row in enumerate(fp):
        if i < selected:
            continue
        print row.rstrip()
        try: 
            selected = selected_idx.next()
        except StopIteration:
            break

    
if __name__ == '__main__':
    
    p = OptionParser(__doc__)

    filter_choices = ("ref", "query", "intersection", "union")
    dialect_choices = ("blast", "coords")

    p.add_option("-f", "--filter", choices=filter_choices,
            dest="filter", default="intersection", 
            help="filter choices: " + str(filter_choices) + " [default: %default]")
    p.add_option("-d", "--dialect", choices=dialect_choices,
            dest="dialect", default=None,
            help="dialect choices: " + str(dialect_choices))

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    blast_file = args[0]

    dialect = opts.dialect
    if not dialect:
        # guess from the suffix
        dialect = "coords" if blast_file.endswith(".coords") else "blast"
        logging.debug("dialect is %s" % dialect)

    filter = opts.filter 

    main(blast_file, filter=filter, dialect=dialect)
