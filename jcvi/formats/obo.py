#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog obo_file

Parses obo_file and plot GO lineage
"""
import sys
import logging

from collections import deque
from functools import partial
from typing import IO, Optional

from goatools.obo_parser import GODag
from jcvi.apps.base import OptionParser

GO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"
SO_URL = (
    "http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/genomic-proteomic/so.obo"
)


def load_GODag(obo_url: str, prt: Optional[IO] = None) -> (GODag, str):
    """
    Load given obo url and returns GODag object.

    Args:
        obo_url (str): URL to the remote OBO file.
        prt (Optional[IO]): IO stream to print verbose information.

    Returns:
        (GODag, str): GODag object that contains the dict, and path to the downloaded OBO file.
    """

    from jcvi.apps.base import download

    so_file = download(obo_url, debug=False)

    return GODag(so_file, prt=prt), so_file


GODag_from_GO = partial(load_GODag, obo_url=GO_URL)
GODag_from_SO = partial(load_GODag, obo_url=SO_URL)


def validate_term(term, so=None, method="verify"):
    """
    Validate an SO term against so.obo
    """
    if so is None:
        so, _ = GODag_from_SO()

    oterm = term
    valid_names = set(x.name for x in so.values())
    if term not in valid_names:
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
        "--term",
        help="Write the parents and children of this query term",
    )

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    (obo_file,) = args

    def description(record):
        level = "level-{:>02}".format(record.level)
        desc = "{} [{}]".format(record.name, record.namespace)
        if record.is_obsolete:
            desc += " obsolete"
        alt_ids = ",".join(record.alt_ids)
        return "\t".join((record.item_id, level, desc, alt_ids))

    g = GODag(obo_file, prt=None)
    header = "\t".join(("#id", "level", "name", "alt_ids"))
    print(header)
    for rec in sorted(set(g.values()), key=lambda x: x.item_id):
        print(description(rec))

    # run a test case
    if opts.term:
        rec = g.query_term(opts.term, verbose=True)
        g.draw_lineage([rec])
