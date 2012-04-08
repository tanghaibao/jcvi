#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for mummerplot. Selecting a subset of queries and references to plot
main features in the dot plot.
"""

import sys
import logging

from optparse import OptionParser

from jcvi.formats.coords import Coords, fromdelta
from jcvi.formats.sizes import Sizes
from jcvi.formats.base import SetFile
from jcvi.apps.base import debug, sh, need_update
debug()


def writeXfile(ids, dict, filename):
    fw = open(filename, "w")
    for q in ids:
        print >> fw, "\t".join(str(x) for x in (q, dict[q], "+"))

    logging.debug("{0} ids written to `{1}`.".format(len(ids), filename))
    fw.close()


def main(args):
    """
    %prog deltafile queryidsfile query.fasta ref.fasta

    Plot one query. Extract the references that have major matches to this
    query. Control "major" by option --querycov.
    """
    p = OptionParser(main.__doc__)
    p.add_option("--querycov", default=.01, type="float",
                 help="Minimum query coverage [default: %default]")
    p.add_option("--all", default=False, action="store_true",
                 help="Plot one pdf file per query [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    deltafile, queryidsfile, queryfasta, reffasta = args
    qsizes = Sizes(queryfasta).mapping
    rsizes = Sizes(reffasta).mapping
    queries = SetFile(queryidsfile)
    querycov = opts.querycov

    if opts.all:
        for q in queries:
            pdffile = plot_some_queries([q], qsizes, rsizes, deltafile, querycov)
            sh("mv {0} {1}.pdf".format(pdffile, q))
    else:
        plot_some_queries(queries, qsizes, rsizes, deltafile, querycov)


def plot_some_queries(queries, qsizes, rsizes, deltafile, querycov):

    Qfile, Rfile = "Qfile", "Rfile"
    coords = Coords(deltafile)
    refs = set()
    for c in coords:
        if c.querycov < querycov:
            continue
        if c.query not in queries:
            continue
        refs.add(c.ref)

    writeXfile(queries, qsizes, Qfile)
    writeXfile(refs, rsizes, Rfile)

    cmd = "mummerplot {0}".format(deltafile)
    cmd += " -Rfile {0} -Qfile {1}".format(Rfile, Qfile)
    cmd += " --postscript --layout"
    sh(cmd)

    prefix = "out"
    cmd = "ps2pdf {0}.ps {0}.pdf".format(prefix)
    sh(cmd)

    return prefix + ".pdf"


if __name__ == '__main__':
    main(sys.argv[1:])
