#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Programatically accessing UniprotKB to get data from a list of queries
"""
import os.path as op
import sys
import time
import logging

from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

from jcvi.formats.base import must_open
from jcvi.apps.base import OptionParser, ActionDispatcher


uniprot_url = "http://www.uniprot.org/uniprot/"

valid_formats = [
    "html",
    "tab",
    "xls",
    "fasta",
    "gff",
    "txt",
    "xml",
    "rdf",
    "list",
    "rss",
]
valid_columns = [
    "citation",
    "clusters",
    "comments",
    "database",
    "domains",
    "domain",
    "ec",
    "id",
    "entry name",
    "existence",
    "families",
    "features",
    "genes",
    "go",
    "go-id",
    "interpro",
    "interactor",
    "keywords",
    "keyword-id",
    "last-modified",
    "length",
    "organism",
    "organism-id",
    "pathway",
    "protein names",
    "reviewed",
    "score",
    "sequence",
    "3d",
    "subcellular locations",
    "taxon",
    "tools",
    "version",
    "virus hosts",
]

valid_column_formats = ["tab", "xls"]
valid_include_formats = ["fasta", "rdf"]


def main():

    actions = (("fetch", "fetch records from uniprot. input is a list of query terms"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fetch(args):
    """
    %prog fetch "query"
        OR
    %prog fetch queries.txt

    Please provide a UniProt compatible `query` to retrieve data. If `query` contains
    spaces, please remember to "quote" it.

    You can also specify a `filename` which contains queries, one per line.

    Follow this syntax <http://www.uniprot.org/help/text-search#text-search-syntax>
    to query any of the documented fields <http://www.uniprot.org/help/query-fields>
    """
    import re
    import csv

    p = OptionParser(fetch.__doc__)

    p.add_option(
        "--format",
        default="tab",
        choices=valid_formats,
        help="download format",
    )
    p.add_option(
        "--columns",
        default="entry name, protein names, genes,organism",
        help="columns to download, if --format is `tab` or `xls`",
    )
    p.add_option(
        "--include",
        default=False,
        action="store_true",
        help="Include isoforms when --format is `fasta` or include `description` when --format is `rdf`.",
    )
    p.add_option(
        "--limit",
        default=10,
        type="int",
        help="Max number of results to retrieve",
    )
    p.add_option(
        "--offset",
        default=0,
        type="int",
        help="Offset of first result, used with --limit",
    )
    p.add_option(
        "--skipcheck",
        default=False,
        action="store_true",
        help="Turn off prompt to check file existence",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (query,) = args
    url_params = {}
    if op.exists(query):
        pf = query.rsplit(".", 1)[0]
        list_of_queries = [row.strip() for row in open(query)]
    else:
        # the query is the search term
        pf = query.strip().strip('"')
        list_of_queries = [pf]
        pf = re.sub(r"\s+", "_", pf)

    assert len(list_of_queries) > 0, "Please provide atleast one input query"

    url_params["format"] = opts.format

    if opts.columns and opts.format in valid_column_formats:
        reader = csv.reader([opts.columns], skipinitialspace=True)
        cols = [col for r in reader for col in r]
        for col in cols:
            assert (
                col in valid_columns
            ), "Column '{0}' is not a valid. Allowed options are {1}".format(
                col, valid_columns
            )
        url_params["columns"] = ",".join(cols)

    if opts.include and opts.format in valid_include_formats:
        url_params["include"] = "yes"

    url_params["limit"] = opts.limit
    url_params["offset"] = opts.offset

    outfile = "{0}.{1}".format(pf, opts.format)

    # If noprompt, will not check file existence
    fw = must_open(outfile, "w", checkexists=True, skipcheck=opts.skipcheck)
    if fw is None:
        return

    seen = set()
    for query in list_of_queries:
        if query in seen:
            logging.error("Duplicate query ({0}) found".format(query))
            continue

        url_params["query"] = query

        data = urlencode(url_params)
        try:
            request = Request(uniprot_url, data)
            response = urlopen(request)
        except (HTTPError, URLError, RuntimeError, KeyError) as e:
            logging.error(e)
            logging.debug("wait 5 seconds to reconnect...")
            time.sleep(5)

        page = response.read()
        if not page:
            logging.error("query `{0}` yielded no results".format(query))
            continue

        print(page, file=fw)

        seen.add(query)

    if seen:
        print(
            "A total of {0} out of {1} queries returned results.".format(
                len(seen), len(list_of_queries)
            ),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
