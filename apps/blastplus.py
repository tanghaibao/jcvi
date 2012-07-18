#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import sys
import logging

from optparse import OptionParser
from subprocess import Popen, PIPE
from multiprocessing import Lock, Pool
from itertools import islice
from Bio import SeqIO

from jcvi.formats.base import must_open
from jcvi.apps.grid import Grid, Jobs
from jcvi.apps.base import ActionDispatcher, debug, set_params, \
        set_grid, set_outfile, sh, mkdir
from jcvi.apps.command import run_formatdb
debug()


blastplus_template = "{0} -query {1} -db {2} -out {3} -outfmt {4} "

# recommanded params for comparing annotations of two genomes, or a self-blast
recommendOptions = " -evalue 1e-5 -num_descriptions 20 -num_alignments 20 "


def blastplus(k, n, bfasta_fn, afasta_fn, out_fh, lock, blast_bin, extra, \
format, grid=False):
    blast_cmd = blastplus_template.\
    format(blast_bin, afasta_fn, bfasta_fn, out_fh.name, format)
    if extra:
        blast_cmd += " " + extra.strip()
    
    if grid: # if run on SGE, only the cmd is needed
        return blast_cmd    
    
    proc = Popen(blast_cmd, stdin=PIPE, stdout=PIPE, shell=True)
    parser = SeqIO.parse(afasta_fn, "fasta")
    for rec in islice(parser, k - 1, None, n):
        SeqIO.write([rec], proc.stdin, "fasta")
    proc.stdin.close()
    
    logging.debug("job <%d> started: %s" % (proc.pid, blast_cmd))
    for row in proc.stdout:
        lock.acquire()
        out_fh.write(row)
        out_fh.flush()
        lock.release()
    logging.debug("job <%d> finished: %s" % (proc.pid, blast_cmd))


def main():
    """
    %prog database.fa query.fa [options]

    Wrapper for NCBI BLAST+.
    """
    p = OptionParser(main.__doc__)

    p.add_option("-a", "-A", dest="cpus", default=1, type="int",
            help="parallelize job to multiple cpus [default: %default]")
    p.add_option("--format", default=" \'6 qseqid sseqid pident length " \
    "mismatch gapopen qstart qend sstart send evalue bitscore\' ", 
            help="0-11, learn more with \"blastp -help\". [default: %default]")
    p.add_option("--path", dest="blast_path", default=None,
            help="specify BLAST+ path including the program name")
    p.add_option("--program", dest="blast_program", default=None,
            help="specify BLAST+ program to use. See complete list here: " \
            "http://www.ncbi.nlm.nih.gov/books/NBK52640/#chapter1.Installation")
    p.add_option("--recommend", default=False, action="store_true",
            help="Use recommended options tuned for comparing " \
            "whole genome annotations [default: %default]")

    set_params(p)
    set_outfile(p)
    set_grid(p)

    opts, args = p.parse_args()

    if len(args) != 2 or opts.blast_program is None:
        sys.exit(p.print_help())

    bfasta_fn, afasta_fn = args
    for fn in (afasta_fn, bfasta_fn):
        assert op.exists(fn)

    afasta_fn = op.abspath(afasta_fn)
    bfasta_fn = op.abspath(bfasta_fn)
    out_fh = must_open(opts.outfile, "w")

    grid = opts.grid
    if grid:
        print >>sys.stderr, "Running jobs on JCVI grid"

    extra = opts.extra
    if opts.recommend:
        extra += recommendOptions

    blast_bin = opts.blast_path or opts.blast_program
    if op.basename(blast_bin)!=opts.blast_program:
        blast_bin = "".join([blast_bin, "/", opts.blast_program])

    blast_program = opts.blast_program
    cpus = opts.cpus
    logging.debug("Dispatch job to %d cpus" % cpus)
    format = opts.format
    
    dbtype = "prot" if op.basename(blast_bin) in ["blastp", "blastx"] \
    else "nucl"
    run_formatdb(infile=bfasta_fn, outfile="t", dbtype=dbtype)

    outdir = "outdir"
    lock = Lock()

    if grid:
        cmds = [blastplus(k + 1, cpus, bfasta_fn, afasta_fn, out_fh, \
                lock, blast_bin, extra, format, grid) for k in xrange(cpus)]
        mkdir(outdir)
        g = Grid(cmds, outfiles=[op.join(outdir, "out.{0}.blastplus").\
                format(i) for i in range(len(cmds))])
        g.run()
        g.writestatus()

    else:
        args = [(k + 1, cpus, bfasta_fn, afasta_fn, out_fh,
                lock, blast_bin, extra, format) for k in xrange(cpus)]
        g = Jobs(target=blastplus, args=args)
        g.run()


if __name__ == '__main__':
    main()
