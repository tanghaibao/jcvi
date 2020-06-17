#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog obo_file

Parses obo_file and plot GO lineage
"""
from __future__ import print_function

import sys

from goatools.obo_parser import GODag
from jcvi.apps.base import OptionParser


if __name__ == "__main__":
    p = OptionParser(__doc__)
    p.add_option(
        "--term", help="Write the parents and children of this query term",
    )

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    (obo_file,) = args

    def description(rec):
        level = "level-{:>02}".format(rec.level)
        description = "{} [{}]".format(rec.name, rec.namespace)
        if rec.is_obsolete:
            description += " obsolete"
        alt_ids = ",".join(rec.alt_ids)
        return "\t".join((rec.item_id, level, description, alt_ids))

    g = GODag(obo_file, prt=None)
    header = "\t".join(("#id", "level", "name", "alt_ids"))
    print(header)
    for rec in sorted(set(g.values()), key=lambda x: x.item_id):
        print(description(rec))

    # run a test case
    if opts.term:
        rec = g.query_term(opts.term, verbose=True)
        g.draw_lineage([rec], verbose=True)
