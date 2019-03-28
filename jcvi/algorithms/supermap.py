#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog infile [options]

This script combines pairwise alignments, sort and filter the alignments.
Infile expect BLAST tabular format (-m8) or nucmer .coords file.

In order to handle dups, we have to run two monotonic chains in both genomes,
first chain using ref, and a second chain using query and we will have options
to keep either the union or the intersection of retain chained alignments from
both genomes, similar to the SUPERMAP algorithm. This operation is symmetrical.
"""
from __future__ import print_function

import sys
import logging

from jcvi.apps.base import OptionParser
from jcvi.formats.blast import BlastLine
from jcvi.formats.coords import CoordsLine
from jcvi.utils.range import Range, range_chain


def BlastOrCoordsLine(filename, filter="ref", dialect="blast", clip=0):
    allowed_filters = ("ref", "query")
    REF, QUERY = range(len(allowed_filters))

    allowed_dialects = ("blast", "coords")
    BLAST, COORDS = range(len(allowed_dialects))

    assert filter in allowed_filters
    filter = allowed_filters.index(filter)

    assert dialect in allowed_dialects
    dialect = allowed_dialects.index(dialect)

    fp = open(filename)
    for i, row in enumerate(fp):
        if row[0] == '#':
            continue
        if dialect == BLAST:
            b = BlastLine(row)
            if filter == QUERY:
                query, start, end = b.query, b.qstart, b.qstop
            else:
                query, start, end = b.subject, b.sstart, b.sstop
        else:
            try:
                b = CoordsLine(row)
            except AssertionError:
                continue

            if filter == QUERY:
                query, start, end = b.query, b.start2, b.end2
            else:
                query, start, end = b.ref, b.start1, b.end1

        if start > end:
            start, end = end, start

        if clip:
            # clip cannot be more than 5% of the range
            r = end - start + 1
            cc = min(.05 * r, clip)
            start = start + cc
            end = end - cc

        yield Range(query, start, end, b.score, i)


def supermap(blast_file, filter="intersection", dialect="blast", clip=0):
    # filter by query
    if filter != "ref":
        logging.debug("filter by query")
        ranges = list(BlastOrCoordsLine(blast_file, filter="query",
            dialect=dialect, clip=clip))

        query_selected, query_score = range_chain(ranges)
        query_idx = set(x.id for x in query_selected)

    # filter by ref
    if filter != "query":
        logging.debug("filter by ref")
        ranges = list(BlastOrCoordsLine(blast_file, filter="ref",
            dialect=dialect, clip=clip))

        ref_selected, ref_score = range_chain(ranges)
        ref_idx = set(x.id for x in ref_selected)

    if filter == "ref":
        selected_idx = ref_idx

    elif filter == "query":
        selected_idx = query_idx

    elif filter == "intersection":
        logging.debug("perform intersection")
        selected_idx = ref_idx & query_idx

    elif filter == "union":
        logging.debug("perform union")
        selected_idx = ref_idx | query_idx

    assert len(selected_idx) != 0

    # selected_idx is in fact the lineno in the BLAST file
    fp = open(blast_file)

    if filter == "intersection":
        tag = ""
    else:
        tag = "." + filter
    supermapfile = blast_file + tag + ".supermap"
    fw = open(supermapfile, "w")

    selected_idx = iter(sorted(selected_idx))
    selected = next(selected_idx)
    for i, row in enumerate(fp):
        if i < selected:
            continue
        print(row.rstrip(), file=fw)
        try:
            selected = next(selected_idx)
        except StopIteration:
            break

    logging.debug("Write output file to `{0}`".format(supermapfile))
    fw.close()

    from jcvi.formats.blast import sort
    ofilter = "ref" if filter == "ref" else "query"
    args = [supermapfile, "--" + ofilter]
    if dialect == "coords":
        args += ["--coords"]

    sort(args)

    return supermapfile


if __name__ == '__main__':

    p = OptionParser(__doc__)

    filter_choices = ("ref", "query", "intersection", "union")
    dialect_choices = ("blast", "coords")
    p.add_option("--filter", choices=filter_choices, default="intersection",
            help="Available filters [default: %default]")
    p.add_option("--dialect", choices=dialect_choices,
            help="Input format [default: guess]")
    p.add_option("--clip", default=0, type="int",
            help="Clip ranges so that to allow minor overlaps [default: %default]")

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    blast_file, = args

    dialect = opts.dialect
    if not dialect:
        # guess from the suffix
        dialect = "coords" if blast_file.endswith(".coords") else "blast"
        logging.debug("dialect is %s" % dialect)

    supermap(blast_file, filter=opts.filter, dialect=dialect, clip=opts.clip)
