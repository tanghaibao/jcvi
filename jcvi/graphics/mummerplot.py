#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for mummerplot. Selecting a subset of queries and references to plot
main features in the dot plot.
"""
import os.path as op
import sys

from ..apps.base import OptionParser, logger, sh
from ..formats.base import SetFile
from ..formats.coords import Coords, filter
from ..formats.sizes import Sizes


def writeXfile(ids, sizes_dict, filename):
    fw = open(filename, "w")
    for q in ids:
        print("\t".join(str(x) for x in (q, sizes_dict[q], "+")), file=fw)

    logger.debug("%d ids written to `%s`.", len(ids), filename)
    fw.close()


def main(args):
    """
    %prog deltafile

    Plot one query. Extract the references that have major matches to this
    query. Control "major" by option --refcov.
    """
    p = OptionParser(main.__doc__)
    p.add_argument("--refids", help="Use subset of contigs in the ref")
    p.add_argument(
        "--refcov",
        default=0.01,
        type=float,
        help="Minimum reference coverage",
    )
    p.add_argument(
        "--all",
        default=False,
        action="store_true",
        help="Plot one pdf file per ref in refidsfile",
    )
    p.add_argument(
        "--color",
        default="similarity",
        choices=("similarity", "direction", "none"),
        help="Color the dots based on",
    )
    p.add_argument(
        "--nolayout",
        default=False,
        action="store_true",
        help="Do not rearrange contigs",
    )
    p.set_align(pctid=0, hitlen=0)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (deltafile,) = args
    reffasta, queryfasta = open(deltafile).readline().split()
    color = opts.color
    layout = not opts.nolayout
    prefix = op.basename(deltafile).split(".")[0]
    qsizes = Sizes(queryfasta).mapping
    rsizes = Sizes(reffasta).mapping

    refs = SetFile(opts.refids) if opts.refids else set(rsizes.keys())
    refcov = opts.refcov
    pctid = opts.pctid
    hitlen = opts.hitlen
    deltafile = filter(
        [deltafile, "--pctid={0}".format(pctid), "--hitlen={0}".format(hitlen)]
    )

    if opts.all:
        for r in refs:
            pdffile = plot_some_queries(
                [r],
                qsizes,
                rsizes,
                deltafile,
                refcov,
                prefix=prefix,
                color=color,
                layout=layout,
            )
            if pdffile:
                sh("mv {0} {1}.pdf".format(pdffile, r))
    else:
        plot_some_queries(
            refs,
            qsizes,
            rsizes,
            deltafile,
            refcov,
            prefix=prefix,
            color=color,
            layout=layout,
        )


def plot_some_queries(
    refs,
    qsizes,
    rsizes,
    deltafile,
    refcov,
    prefix="out",
    color="similarity",
    layout=True,
):

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
        logger.debug("Empty - %s vs. %s", queries, refs)
        return None

    if not layout:
        queries = sorted(queries)
        refs = sorted(refs)

    writeXfile(queries, qsizes, Qfile)
    writeXfile(refs, rsizes, Rfile)

    cmd = "mummerplot {0}".format(deltafile)
    cmd += " -Rfile {0} -Qfile {1}".format(Rfile, Qfile)
    cmd += " --postscript -p {0}".format(prefix)
    if layout:
        cmd += " --layout"
    if color == "similarity":
        cmd += " --color"
    elif color == "none":
        cmd += " --nocolor"
    sh(cmd)

    cmd = "ps2pdf {0}.ps {0}.pdf".format(prefix)
    sh(cmd)

    return prefix + ".pdf"


if __name__ == "__main__":
    main(sys.argv[1:])
