#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog obo_file

Parses obo_file and plot GO lineage
"""
from __future__ import print_function

import sys
import logging

from collections import deque

from goatools.obo_parser import GODag
from jcvi.apps.base import OptionParser


def load_GODag():
    """
    OBO file retrieved from http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/genomic-proteomic/so.obo
    """
    from jcvi.apps.base import download

    so_file_url = "http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/genomic-proteomic/so.obo"
    so_file = download(so_file_url, debug=False)

    return GODag(so_file)


def validate_term(term, so=None, method="verify"):
    """
    Validate an SO term against so.obo
    """
    if so is None:
        so = load_GODag()

    oterm = term
    valid_names = set(x.name for x in so.values())
    if term not in so.valid_names:
        if "resolve" in method:
            if "_" in term:
                tparts = deque(term.split("_"))
                tparts.pop() if "prefix" in method else tparts.popleft()
                nterm = "_".join(tparts).strip()
                term = validate_term(nterm, so=so, method=method)
            if term is None:
                return None
        else:
            logging.error("Term `{0}` does not exist".format(term))
            sys.exit(1)

    if oterm != term:
        logging.debug("Resolved term `{0}` to `{1}`".format(oterm, term))
    return term


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
