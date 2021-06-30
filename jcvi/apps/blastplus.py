#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import sys
import logging

from multiprocessing import Lock

from jcvi.formats.base import must_open, split
from jcvi.apps.grid import Jobs
from jcvi.apps.align import run_formatdb
from jcvi.apps.base import OptionParser, Popen


def blastplus(out_fh, cmd, query, lock):
    cmd += " -query {0}".format(query)
    proc = Popen(cmd)

    logging.debug("job <%d> started: %s" % (proc.pid, cmd))
    for row in proc.stdout:
        if row[0] == "#":
            continue
        lock.acquire()
        out_fh.write(row)
        out_fh.flush()
        lock.release()
    logging.debug("job <%d> finished" % proc.pid)


def main():
    """
    %prog database.fa query.fa [options]

    Wrapper for NCBI BLAST+.
    """
    p = OptionParser(main.__doc__)

    p.add_option(
        "--format",
        default=" '6 qseqid sseqid pident length "
        "mismatch gapopen qstart qend sstart send evalue bitscore' ",
        help='0-11, learn more with "blastp -help"',
    )
    p.add_option(
        "--path",
        dest="blast_path",
        default=None,
        help="specify BLAST+ path including the program name",
    )
    p.add_option(
        "--prog",
        dest="blast_program",
        default="blastp",
        help="specify BLAST+ program to use. See complete list here: "
        "http://www.ncbi.nlm.nih.gov/books/NBK52640/#chapter1.Installation",
    )
    p.set_align(evalue=0.01)
    p.add_option(
        "--best",
        default=1,
        type="int",
        help="Only look for best N hits",
    )
    p.set_cpus()
    p.add_option(
        "--nprocs",
        default=1,
        type="int",
        help="number of BLAST processes to run in parallel. "
        + "split query.fa into `nprocs` chunks, "
        + "each chunk uses -num_threads=`cpus`",
    )
    p.set_params()
    p.set_outfile()
    opts, args = p.parse_args()

    if len(args) != 2 or opts.blast_program is None:
        sys.exit(not p.print_help())

    bfasta_fn, afasta_fn = args
    for fn in (afasta_fn, bfasta_fn):
        assert op.exists(fn)

    afasta_fn = op.abspath(afasta_fn)
    bfasta_fn = op.abspath(bfasta_fn)
    out_fh = must_open(opts.outfile, "w")

    extra = opts.extra
    blast_path = opts.blast_path
    blast_program = opts.blast_program

    blast_bin = blast_path or blast_program
    if op.basename(blast_bin) != blast_program:
        blast_bin = op.join(blast_bin, blast_program)

    nprocs, cpus = opts.nprocs, opts.cpus
    if nprocs > 1:
        logging.debug("Dispatch job to %d processes" % nprocs)
        outdir = "outdir"
        fs = split([afasta_fn, outdir, str(nprocs)])
        queries = fs.names
    else:
        queries = [afasta_fn]

    dbtype = "prot" if op.basename(blast_bin) in ("blastp", "blastx") else "nucl"

    db = bfasta_fn
    if dbtype == "prot":
        nin = db + ".pin"
    else:
        nin00 = db + ".00.nin"
        nin = nin00 if op.exists(nin00) else (db + ".nin")

    run_formatdb(infile=db, outfile=nin, dbtype=dbtype)

    lock = Lock()

    blastplus_template = "{0} -db {1} -outfmt {2}"
    blast_cmd = blastplus_template.format(blast_bin, bfasta_fn, opts.format)
    blast_cmd += " -evalue {0} -max_target_seqs {1}".format(opts.evalue, opts.best)
    blast_cmd += " -num_threads {0}".format(cpus)
    if extra:
        blast_cmd += " " + extra.strip()

    args = [(out_fh, blast_cmd, query, lock) for query in queries]
    g = Jobs(target=blastplus, args=args)
    g.run()


if __name__ == "__main__":
    main()
