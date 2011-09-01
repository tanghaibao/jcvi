#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import sys
import math
import logging

from optparse import OptionParser
from subprocess import Popen, PIPE
from multiprocessing import Process, Lock

from jcvi.apps.grid import Grid
from jcvi.apps.base import ActionDispatcher, debug, set_params, set_grid
debug()


blast_fields = "query,subject,pctid,hitlen,nmismatch,ngaps,"\
        "qstart,qstop,sstart,sstop,evalue,score"

lastz_fields = "name2,name1,identity,nmismatch,ngap,"\
        "start2+,end2+,strand2,start1,end1,strand1,score"

# conversion between blastz and ncbi is taken from Kent src
# src/lib/blastOut.c
# this is not rigorous definition of e-value (assumes human genome) !!
blastz_score_to_ncbi_bits = lambda bz_score: bz_score * 0.0205


def blastz_score_to_ncbi_expectation(bz_score):
    bits = blastz_score_to_ncbi_bits(bz_score)
    log_prob = -bits * 0.693147181
    # this number looks like.. human genome?
    return 3.0e9 * math.exp(log_prob)


def lastz_to_blast(row):
    """
    Convert the lastz tabular to the blast tabular, see headers above
    Obsolete after LASTZ version 1.02.40
    """
    atoms = row.strip().split("\t")
    name1, name2, coverage, identity, nmismatch, ngap, \
            start1, end1, strand1, start2, end2, strand2, score = atoms
    identity = identity.replace("%", "")
    hitlen = coverage.split("/")[1]
    score = float(score)
    same_strand = (strand1 == strand2)
    if not same_strand:
        start2, end2 = end2, start2

    evalue = blastz_score_to_ncbi_expectation(score)
    score = blastz_score_to_ncbi_bits(score)
    evalue, score = "%.2g" % evalue, "%.1f" % score
    return "\t".join((name1, name2, identity, hitlen, nmismatch, ngap, \
            start1, end1, start2, end2, evalue, score))


def lastz(k, n, bfasta_fn, afasta_fn, out_fh, lock, lastz_path, extra,
        blastline=True, mask=False, grid=False):
    lastz_bin = lastz_path or "lastz"

    ref_tags = ["multiple", "nameparse=darkspace"]
    qry_tags = ["nameparse=darkspace", "subsample=%d/%d" % (k, n)]
    if not mask:
        ref_tags.append("unmask")
        qry_tags.append("unmask")

    ref_tags = ",".join(ref_tags)
    qry_tags = ",".join(qry_tags)

    lastz_cmd = "%s --ambiguous=iupac %s[%s] %s[%s] %s"
    lastz_cmd %= (lastz_bin, bfasta_fn, ref_tags, afasta_fn, qry_tags, extra)

    if blastline:
        #lastz_cmd += " --format=general-:%s" % lastz_fields
        # The above conversion is no longer necessary after LASTZ v1.02.40
        # (of which I contributed a patch)
        lastz_cmd += " --format=BLASTN-"

    if grid:  # if run on SGE, only the cmd is needed
        return lastz_cmd

    proc = Popen(lastz_cmd, bufsize=1, stdout=PIPE, shell=True)

    logging.debug("job <%d> started: %s" % (proc.pid, lastz_cmd))
    for row in proc.stdout:
        #if blastline: row = lastz_to_blast(row)
        lock.acquire()
        out_fh.write(row)
        out_fh.flush()
        lock.release()
    logging.debug("job <%d> finished" % proc.pid)


def main():
    """
    %prog database.fa query.fa [options]

    Run LASTZ similar to the BLAST interface, and generates -m8 tabular format
    """
    p = OptionParser(main.__doc__)

    p.add_option("-o", dest="outfile",
            help="output [default: stdout]")
    p.add_option("-m", dest="blastline", default=True, action="store_false",
            help="don't generate BLAST tabular format (as -m8) [default: m8]")
    p.add_option("-a", "-A", dest="cpus", default=1, type="int",
            help="parallelize job to multiple cpus [default: %default]")
    p.add_option("--path", dest="lastz_path", default=None,
            help="specify LASTZ path")
    p.add_option("--mask", dest="mask", default=False, action="store_true",
            help="treat lower-case letters as mask info [default: %default]")
    set_params(p)
    set_grid(p)

    opts, args = p.parse_args()

    if len(args) != 2:
        sys.exit(p.print_help())

    bfasta_fn, afasta_fn = args
    for fn in (afasta_fn, bfasta_fn):
        assert op.exists(fn)

    afasta_fn = op.abspath(afasta_fn)
    bfasta_fn = op.abspath(bfasta_fn)
    out_fh = file(opts.outfile, "w") if opts.outfile else sys.stdout

    grid = opts.grid
    if grid:
        print >>sys.stderr, "Running jobs on JCVI grid"

    extra = opts.extra

    lastz_path = opts.lastz_path
    cpus = opts.cpus
    logging.debug("Dispatch job to %d cpus" % cpus)

    lock = Lock()

    if grid:
        cmds = []
        for k in xrange(cpus):
            lastz_cmd = lastz(k + 1, cpus, bfasta_fn, afasta_fn, out_fh,
                    lock, lastz_path, extra, blastline=opts.blastline,
                    mask=opts.mask, grid=grid)
            cmds.append(lastz_cmd)

        g = Grid(cmds, outfiles=["lastz.out.{0}".\
                format(i) for i in range(len(cmds))])
        g.run()
        g.writestatus()

    else:
        processes = []
        for k in xrange(cpus):
            pi = Process(target=lastz, args=(k + 1, cpus,
                bfasta_fn, afasta_fn, out_fh,
                lock, lastz_path, extra, opts.blastline, opts.mask))
            pi.start()
            processes.append(pi)

        for pi in processes:
            pi.join()


if __name__ == '__main__':
    main()
