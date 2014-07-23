#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for LAST program.

<http://last.cbrc.jp>
"""

import os.path as op
import sys
import logging

from itertools import islice
from multiprocessing import Lock
from Bio import SeqIO

from jcvi.utils.cbook import depends
from jcvi.apps.grid import Jobs
from jcvi.formats.base import must_open
from jcvi.apps.base import OptionParser, sh, Popen, PIPE


@depends
def run_lastdb(infile=None, outfile=None, mask=False, lastdb_bin="lastdb"):
    outfilebase = outfile.rsplit(".", 1)[0]
    mask = "-c " if mask else ""
    cmd = "{0} {1}{2} {3}".format(lastdb_bin, mask, outfilebase, infile)
    sh(cmd)


def last(k, n, out_fh, cmd, query, lock):

    proc = Popen(cmd, stdin=PIPE)

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
    p.add_option("--path", help="specify LAST path")
    p.add_option("--mask", default=False, action="store_true",
                 help="invoke -c in lastdb [default: %default]")
    p.add_option("--format", default="blast", choices=supported_formats,
                 help="Output format [default: %default]")
    p.add_option("--eval", default=False, action="store_true",
                 help="Use lastex to recalculate E-value [default: %default]")
    p.set_cpus(cpus=32)
    p.set_params()
    p.set_outfile()

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
    run_lastdb(infile=subject, outfile=subjectdb + ".prj", mask=opts.mask, \
              lastdb_bin=lastdb_bin)

    cpus = opts.cpus
    logging.debug("Dispatch job to {0} cpus".format(cpus))

    oappend = False
    if opts.format == "maf":
        cmd = 'echo "##maf version=1"'
        sh(cmd, outfile=opts.outfile)
        oappend = True

    u = 2 if opts.mask else 0
    cmd = "{0} -u {1}".format(lastal_bin, u)
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

    out_fh = must_open(opts.outfile, "w", checkexists=True, oappend=oappend)

    if out_fh is None:
        return

    lock = Lock()
    args = [(k + 1, cpus, out_fh, cmd, query, lock) \
                    for k in xrange(cpus)]
    g = Jobs(target=last, args=args)
    g.run()


if __name__ == '__main__':
    main(sys.argv[1:])
