#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for LAST program.

<http://last.cbrc.jp>
"""

import sys

from optparse import OptionParser

from jcvi.utils.cbook import depends
from jcvi.apps.base import ActionDispatcher, debug, sh, set_outfile
debug()


def main():

    actions = (
        ('last', 'run LAST by calling LASTDB, LASTAL and LASTEX'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


@depends
def run_lastdb(infile=None, outfile=None):
    outfilebase = outfile.rsplit(".", 1)[0]
    cmd = "lastdb -c {0} {1}".format(outfilebase, infile)
    sh(cmd)


def last(args):
    """
    %prog last database.fasta query.fasta


    Run LAST by calling LASTDB, LASTAL and LASTEX.
    """

    supported_formats = ("tab", "maf")

    p = OptionParser(last.__doc__)
    p.add_option("--format", default="tab", choices=supported_formats,
                 help="Output format, one of {0} [default: %default]".\
                      format("|".join(supported_formats)))
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    subject, query = args
    subjectdb = subject.rsplit(".", 1)[0]
    querydb = query.rsplit(".", 1)[0]

    run_lastdb(infile=subject, outfile=subjectdb + ".prj")
    run_lastdb(infile=query, outfile=querydb + ".prj")

    cmd = "lastal"
    if opts.format == "tab":
        cmd += " -f 0"
    cmd += " {0} {1}".format(subjectdb, query)
    cmd += " | lastex {0}.prj {1}.prj -".format(subjectdb, querydb)
    sh(cmd)


if __name__ == '__main__':
    main()
