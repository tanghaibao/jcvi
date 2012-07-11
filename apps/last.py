#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for LAST program.

<http://last.cbrc.jp>
"""

import os.path as op
import sys
import logging

from optparse import OptionParser
from subprocess import Popen, PIPE

from Bio import SeqIO
from itertools import islice
from multiprocessing import Lock

from jcvi.utils.cbook import depends
from jcvi.apps.grid import Jobs
from jcvi.formats.base import must_open
from jcvi.apps.base import debug, sh, set_outfile, set_params
debug()


@depends
def run_lastdb(infile=None, outfile=None, lastdb_bin="lastdb"):
    outfilebase = outfile.rsplit(".", 1)[0]
    cmd = "{0} {1} {2}".format(lastdb_bin, outfilebase, infile)
    sh(cmd)


def last(k, n, out_fh, cmd, query, lock):

    proc = Popen(cmd, stdin=PIPE, stdout=PIPE, shell=True)

    parser = SeqIO.parse(query, "fasta")
    for rec in islice(parser, k - 1, None, n):
        SeqIO.write([rec], proc.stdin, "fasta")
    proc.stdin.close()

    logging.debug("job <%d> started: %s" % (proc.pid, cmd))
    for row in proc.stdout:
        if row[0] == '#':
            continue
        lock.acquire()
        out_fh.write(row)
        out_fh.flush()
        lock.release()
    logging.debug("job <%d> finished" % proc.pid)


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

    set_params(p)
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    subject, query = args
    if opts.eval and opts.cpus > 1:
        raise Exception, "Option --eval cannnot work with multiple threads"

    path = opts.path
    getpath = lambda x: op.join(path, x) if path else x
    lastdb_bin = getpath("lastdb")
    lastal_bin = getpath("lastal")
    lastex_bin = getpath("lastex")

    subjectdb = subject.rsplit(".", 1)[0]
    run_lastdb(infile=subject, outfile=subjectdb + ".prj", lastdb_bin=lastdb_bin)

    cpus = opts.cpus
    logging.debug("Dispatch job to {0} cpus".format(cpus))

    if opts.format == "maf":
        cmd = 'echo "##maf version=1"'
        sh(cmd)

    cmd = "{0} -u 0".format(lastal_bin)
    f = supported_formats.index(opts.format)
    cmd += " -f {0}".format(f)
    cmd += " {0} -".format(subjectdb)

    extra = opts.extra
    if extra:
        cmd += " " + extra

    if opts.eval:
        querydb = query.rsplit(".", 1)[0]
        run_lastdb(infile=query, outfile=querydb + ".prj")

        cmd += " | {0} {1}.prj {2}.prj -".format(lastex_bin, subjectdb, querydb)

    out_fh = must_open(opts.outfile, "w")
    lock = Lock()

    args = [(k + 1, cpus, out_fh, cmd, query, lock) \
                    for k in xrange(cpus)]
    g = Jobs(target=last, args=args)
    g.run()


if __name__ == '__main__':
    main(sys.argv[1:])
