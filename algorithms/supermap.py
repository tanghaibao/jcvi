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
from jcvi.utils.range import Range, range_chain
from jcvi.apps.base import debug
debug()


def main(blast_file, opts):
    # filter by query
    if opts.filter!="ref":
        logging.debug("now filter by query")
        ranges = []
        fp = open(blast_file)

        for i, row in enumerate(fp):
            b = BlastLine(row)
            ranges.append(Range(b.query, b.qstart, b.qstop, b.score, i))
        query_selected, query_score = range_chain(ranges)
        query_idx = set(x.id for x in query_selected)

    # filter by ref
    if opts.filter!="query":
        logging.debug("now filter by ref")
        ranges = []
        fp = open(blast_file)

        for i, row in enumerate(fp):
            b = BlastLine(row)
            ranges.append(Range(b.subject, b.sstart, b.sstop, b.score, i))
        ref_selected, ref_score = range_chain(ranges)
        ref_idx = set(x.id for x in ref_selected)

    if opts.filter=="ref":
        selected_idx = ref_idx

    elif opts.filter=="query":
        selected_idx = query_idx

    elif opts.filter=="intersection":
        logging.debug("perform intersection")
        selected_idx = ref_idx & query_idx

    elif opts.filter=="union":
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
        print row.strip()
        try: 
            selected = selected_idx.next()
        except StopIteration:
            break

    
if __name__ == '__main__':
    
    p = OptionParser(__doc__)

    filter_choices = ("ref", "query", "intersection", "union")
    p.add_option("-f", "--filter", choices=filter_choices,
            dest="filter", default="intersection", 
            help="filter choices: " + str(filter_choices) + " [default: %default]")

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    blast_file = args[0]

    main(blast_file, opts)
