#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run through NCBI vecscreen on a local machine.
"""

import sys

from optparse import OptionParser

from jcvi.utils.cbook import depends 
from jcvi.utils.range import range_merge
from jcvi.formats.blast import BlastLine
from jcvi.apps.base import ActionDispatcher, debug, download, sh
debug()


def main():

    actions = (
        ('blast', 'run BLASTN against UniVec'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


@depends
def run_formatdb(infile=None, outfile=None):
    cmd = "formatdb -i {0} -p F".format(infile)
    sh(cmd)


@depends
def run_blastall(infile=None, outfile=None):
    """
    BLASTN parameters reference:
    http://www.ncbi.nlm.nih.gov/VecScreen/VecScreen_docs.html
    """
    cmd = 'blastall -p blastn -i {0}'.format(infile) 
    cmd += ' -d UniVec_Core -q -5 -G 3 -E 3 -F "m D"'
    cmd += ' -e 700 -Y 1.75e12 -m 8 -o {0}'.format(outfile)
    sh(cmd)


@depends
def run_supermap(infile=None, outfile=None):
    from jcvi.algorithms.supermap import supermap

    supermap(fastablast, filter="query")


def blast(args):
    """
    %prog blast fastafile 

    Run BLASTN against UniVec.
    """
    p = OptionParser(blast.__doc__)
    p.add_option("--dist", dest="dist", default=100,
            help="Merge adjacent HSPs separated by [default: %default]")
    opts, args = p.parse_args(args)

    fastafile, = args
    univec = download("ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core")
    univecdb = univec + ".nin"
    run_formatdb(infile=univec, outfile=univecdb)
    fastablast = fastafile + ".uniblast"
    run_blastall(infile=fastafile, outfile=fastablast)
    supermapblast = fastablast + ".query.supermap"
    run_supermap(infile=fastablast, outfile=supermapblast)

    fp = open(supermapblast)
    ranges = []
    for row in fp:
        b = BlastLine(row)
        ranges.append((b.query, b.qstart, b.qstop))

    merged_ranges = range_merge(ranges, dist=opts.dist)
    for seqid, start, end in merged_ranges:
        print "\t".join(str(x) for x in (seqid, start, end, "vector"))


if __name__ == '__main__':
    main()
