#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for LAST program.

<http://last.cbrc.jp>
"""

import sys

from optparse import OptionParser

from jcvi.utils.cbook import depends
from jcvi.apps.base import debug, sh, set_outfile
debug()


@depends
def run_lastdb(infile=None, outfile=None):
    outfilebase = outfile.rsplit(".", 1)[0]
    cmd = "lastdb -c {0} {1}".format(outfilebase, infile)
    sh(cmd)


def main(args):
    """
    %prog database.fasta query.fasta


    Run LAST by calling LASTDB, LASTAL and LASTEX.
    """

    supported_formats = ("tab", "maf", "blast")

    p = OptionParser(main.__doc__)
    p.add_option("-a", "-A", dest="cpus", default=1, type="int",
            help="parallelize job to multiple cpus [default: %default]")
    p.add_option("--path", help="specify LAST path")
    p.add_option("--format", default="blast", choices=supported_formats,
                 help="Output format, one of {0} [default: %default]".\
                      format("|".join(supported_formats)))
    p.add_option("--eval", default=False, action="store_true",
                 help="Use lastex to recalculate E-value [default: %default]")
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    subject, query = args
    if opts.eval and opts.cpus > 1:
        raise Exception, "Option --eval cannnot work with multiple threads"

    subjectdb = subject.rsplit(".", 1)[0]
    querydb = query.rsplit(".", 1)[0]

    run_lastdb(infile=subject, outfile=subjectdb + ".prj")

    if opts.format == "maf":
        cmd = 'echo "##maf version=1"'
        sh(cmd)

    cmd = "lastal -u 0"
    f = supported_formats.index(opts.format)
    cmd += " -f {0}".format(f)
    cmd += " {0} {1}".format(subjectdb, query)

    if opts.eval:
        run_lastdb(infile=query, outfile=querydb + ".prj")
        cmd += " | lastex {0}.prj {1}.prj -".format(subjectdb, querydb)

    sh(cmd)


if __name__ == '__main__':
    main(sys.argv[1:])
