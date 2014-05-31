#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for mummerplot. Selecting a subset of queries and references to plot
main features in the dot plot.
"""

import sys
import logging

from jcvi.formats.coords import Coords, filter
from jcvi.formats.sizes import Sizes
from jcvi.formats.base import SetFile
from jcvi.apps.base import OptionParser, sh


def writeXfile(ids, dict, filename):
    fw = open(filename, "w")
    for q in ids:
        print >> fw, "\t".join(str(x) for x in (q, dict[q], "+"))

    logging.debug("{0} ids written to `{1}`.".format(len(ids), filename))
    fw.close()


def main(args):
    """
    %prog deltafile refidsfile query.fasta ref.fasta

    Plot one query. Extract the references that have major matches to this
    query. Control "major" by option --refcov.
    """
    p = OptionParser(main.__doc__)
    p.add_option("--refcov", default=.01, type="float",
                 help="Minimum reference coverage [default: %default]")
    p.add_option("--all", default=False, action="store_true",
                 help="Plot one pdf file per ref in refidsfile [default: %default]")
    p.set_align(pctid=96, hitlen=500)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    deltafile, refidsfile, queryfasta, reffasta = args
    qsizes = Sizes(queryfasta).mapping
    rsizes = Sizes(reffasta).mapping
    refs = SetFile(refidsfile)
    refcov = opts.refcov
    pctid = opts.pctid
    hitlen = opts.hitlen
    deltafile = filter([deltafile, "--pctid={0}".format(pctid),
                        "--hitlen={0}".format(hitlen)])

    if opts.all:
        for r in refs:
            pdffile = plot_some_queries([r], qsizes, rsizes, deltafile, refcov)
            if pdffile:
                sh("mv {0} {1}.pdf".format(pdffile, r))
    else:
        plot_some_queries(refs, qsizes, rsizes, deltafile, refcov)


def plot_some_queries(refs, qsizes, rsizes, deltafile, refcov):

    Qfile, Rfile = "Qfile", "Rfile"
    coords = Coords(deltafile)
    queries = set()
    for c in coords:
        if c.refcov < refcov:
            continue
        if c.ref not in refs:
            continue
        queries.add(c.query)

    if not queries or not refs:
        logging.debug("Empty - {0} vs. {1}".format(queries, refs))
        return None

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
