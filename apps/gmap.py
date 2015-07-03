#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run GMAP/GSNAP commands. GMAP/GSNAP manual:

<http://research-pub.gene.com/gmap/src/README>
"""

import os.path as op
import sys
import logging

from jcvi.formats.sam import get_prefix
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh, \
            get_abs_path, backup


def main():

    actions = (
        ('index', 'wraps gmap_build'),
        ('align', 'wraps gsnap'),
        ('gmap', 'wraps gmap'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def check_index(dbfile, go=True):
    dbfile = get_abs_path(dbfile)
    dbdir, filename = op.split(dbfile)
    if not dbdir:
        dbdir = "."
    dbname = filename.rsplit(".", 1)[0]
    safile = op.join(dbdir, "{0}/{0}.genomecomp".format(dbname))
    if dbname == filename:
        dbname = filename + ".db"
    if not go:
        return dbdir, dbname

    if need_update(dbfile, safile):
        cmd = "gmap_build -D {0} -d {1} {2}".format(dbdir, dbname, filename)
        sh(cmd)
    else:
        logging.error("`{0}` exists. `gmap_build` already run.".format(safile))

    return dbdir, dbname


def index(args):
    """
    %prog index database.fasta
`
    Wrapper for `gmap_build`. Same interface.
    """
    p = OptionParser(index.__doc__)
    p.add_option("--supercat", default=False, action="store_true",
                 help="Concatenate reference to speed up alignment")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    dbfile, = args
    supercat = opts.supercat
    pf = dbfile.rsplit(".", 1)[0]
    if supercat:
        supercatfile = pf + ".supercat"
        if need_update(dbfile, supercatfile):
            cmd = "tGBS-Generate_Pseudo_Genome.pl"
            cmd += " -f {0} -o {1}".format(dbfile, supercatfile)
            sh(cmd)
        dbfile = supercatfile + ".fasta"
        coordsfile = supercatfile + ".coords"
        # Rename .coords file since gmap_build will overwrite it
        coordsbak = backup(coordsfile)

    check_index(dbfile)
    if supercat:
        sh("mv {0} {1}".format(coordsbak, coordsfile))


def gmap(args):
    """
    %prog gmap database.fasta fastafile

    Wrapper for `gmap`.
    """
    p = OptionParser(gmap.__doc__)
    p.add_option("--cross", default=False, action="store_true",
                 help="Cross-species alignment")
    p.add_option("--npaths", default=0, type="int",
                 help="Maximum number of paths to show."
                 " If set to 0, prints two paths if chimera"
                 " detected, else one.")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    dbfile, fastafile = args
    assert op.exists(dbfile) and op.exists(fastafile)
    prefix = get_prefix(fastafile, dbfile)
    logfile = prefix + ".log"
    gmapfile = prefix + ".gmap.gff3"

    if not need_update((dbfile, fastafile), gmapfile):
        logging.error("`{0}` exists. `gmap` already run.".format(gmapfile))
    else:
        dbdir, dbname = check_index(dbfile)
        cmd = "gmap -D {0} -d {1}".format(dbdir, dbname)
        cmd += " -f 2 --intronlength=100000"  # Output format 2
        cmd += " -t {0}".format(opts.cpus)
        cmd += " --npaths {0}".format(opts.npaths)
        if opts.cross:
            cmd += " --cross-species"
        cmd += " " + fastafile

        sh(cmd, outfile=gmapfile, errfile=logfile)

    return gmapfile, logfile


def align(args):
    """
    %prog align database.fasta read1.fq read2.fq

    Wrapper for `gsnap` single-end or paired-end, depending on the number of
    args.
    """
    from jcvi.formats.fastq import guessoffset

    p = OptionParser(align.__doc__)
    p.add_option("--rnaseq", default=False, action="store_true",
                 help="Input is RNA-seq reads, turn splicing on")
    p.add_option("--snp", default=False, action="store_true",
                 help="Call SNPs after GSNAP")
    p.add_option("--outdir", default=".",
                 help="Output alignment files in directory")
    p.set_home("eddyyeh")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) == 2:
        logging.debug("Single-end alignment")
    elif len(args) == 3:
        logging.debug("Paired-end alignment")
    else:
        sys.exit(not p.print_help())

    dbfile, readfile = args[:2]
    outdir = opts.outdir
    assert op.exists(dbfile) and op.exists(readfile)
    prefix = get_prefix(readfile, dbfile)
    logfile = op.join(outdir, prefix + ".log")
    gsnapfile = op.join(outdir, prefix + ".gsnap")
    nativefile = gsnapfile.rsplit(".", 1)[0] + ".unique.native"

    if not need_update((dbfile, readfile), gsnapfile):
        logging.error("`{0}` exists. `gsnap` already run.".format(gsnapfile))
    else:
        dbdir, dbname = check_index(dbfile)
        cmd = "gsnap -D {0} -d {1}".format(dbdir, dbname)
        cmd += " -B 5 -m 0.1 -i 2 -n 3"  # memory, mismatch, indel penalty, nhits
        if opts.rnaseq:
            cmd += " -N 1"
        cmd += " -t {0}".format(opts.cpus)
        cmd += " --gmap-mode none --nofails"
        if readfile.endswith(".gz"):
            cmd += " --gunzip"
        try:
            offset = "sanger" if guessoffset([readfile]) == 33 else "illumina"
            cmd += " --quality-protocol {0}".format(offset)
        except AssertionError:
            pass
        cmd += " " + " ".join(args[1:])
        sh(cmd, outfile=gsnapfile, errfile=logfile)

    if opts.snp:
        EYHOME = opts.eddyyeh_home
        if need_update(gsnapfile, nativefile):
            cmd = op.join(EYHOME, "convert2native.pl")
            cmd += " --gsnap {0} -o {1}".format(gsnapfile, nativefile)
            cmd += " -proc {0}".format(opts.cpus)
            sh(cmd)

    return gsnapfile, logfile


if __name__ == '__main__':
    main()
