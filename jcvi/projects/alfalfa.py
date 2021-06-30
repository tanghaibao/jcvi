#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Random collection of scripts associated with alfalfa assembly.
"""

import sys

from jcvi.formats.bed import Bed, fastaFromBed
from jcvi.graphics.mummerplot import main as mummerplot_main
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


def main():

    actions = (("nucmer", "select specific chromosome region based on MTR mapping"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def nucmer(args):
    """
    %prog nucmer mappings.bed MTR.fasta assembly.fasta chr1 3

    Select specific chromosome region based on MTR mapping. The above command
    will extract chr1:2,000,001-3,000,000.
    """
    p = OptionParser(nucmer.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 5:
        sys.exit(not p.print_help())

    mapbed, mtrfasta, asmfasta, chr, idx = args
    idx = int(idx)
    m1 = 1000000
    bedfile = "sample.bed"
    bed = Bed()
    bed.add("\t".join(str(x) for x in (chr, (idx - 1) * m1, idx * m1)))
    bed.print_to_file(bedfile)

    cmd = "intersectBed -a {0} -b {1} -nonamecheck -sorted | cut -f4".format(
        mapbed, bedfile
    )
    idsfile = "query.ids"
    sh(cmd, outfile=idsfile)

    sfasta = fastaFromBed(bedfile, mtrfasta)
    qfasta = "query.fasta"
    cmd = "faSomeRecords {0} {1} {2}".format(asmfasta, idsfile, qfasta)
    sh(cmd)

    cmd = "nucmer {0} {1}".format(sfasta, qfasta)
    sh(cmd)

    mummerplot_main(["out.delta", "--refcov=0"])
    sh("mv out.pdf {0}.{1}.pdf".format(chr, idx))


if __name__ == "__main__":
    main()
