#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Finishing pipeline, starting with a phase1/2 BAC. The pipeline ideally should
include the following components

+ BLAST against the Illumina contigs to fish out additional seqs
+ Use minimus2 to combine the contigs through overlaps
+ Map the mates to the contigs and perform scaffolding
+ Base corrections using Illumina reads
"""

import os
import os.path as op
import sys
import logging

from optparse import OptionParser

from jcvi.formats.contig import ContigFile
from jcvi.formats.fasta import Fasta, SeqIO, gaps, format
from jcvi.utils.cbook import depends
from jcvi.assembly.base import n50
from jcvi.apps.command import run_megablast
from jcvi.apps.base import ActionDispatcher, debug, sh, mkdir, is_newer_file
debug()


def main():

    actions = (
        ('overlap', 'build larger contig set by fishing additional seqs'),
        ('overlapbatch', 'call overlap on a set of sequences'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


@depends
def run_gapsplit(infile=None, outfile=None):
    gaps([infile, "--split"])
    return outfile


def overlapbatch(args):
    """
    %prog overlapbatch ctgfasta poolfasta

    Fish out the sequences in `poolfasta` that overlap with `ctgfasta`.
    Mix and combine using `minimus2`.
    """
    p = OptionParser(overlap.__doc__)
    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())

    ctgfasta, poolfasta = args
    f = Fasta(ctgfasta)
    for k, rec in f.iteritems_ordered():
        fastafile = k + ".fasta"
        fw = open(fastafile, "w")
        SeqIO.write([rec], fw, "fasta")
        fw.close()

        overlap([fastafile, poolfasta])


def overlap(args):
    """
    %prog overlap ctgfasta poolfasta

    Fish out the sequences in `poolfasta` that overlap with `ctgfasta`.
    Mix and combine using `minimus2`.
    """
    p = OptionParser(overlap.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ctgfasta, poolfasta = args
    prefix = ctgfasta.split(".")[0]
    rid = list(Fasta(ctgfasta).iterkeys())
    assert len(rid) == 1, "Use overlapbatch() to improve multi-FASTA file"

    rid = rid[0]
    splitctgfasta = ctgfasta.rsplit(".", 1)[0] + ".split.fasta"
    ctgfasta = run_gapsplit(infile=ctgfasta, outfile=splitctgfasta)

    # Run BLAST
    blastfile = ctgfasta + ".blast"
    run_megablast(infile=ctgfasta, outfile=blastfile, db=poolfasta)

    # Extract contigs and merge using minimus2
    closuredir = prefix + ".closure"
    closure = False
    if not op.exists(closuredir) or is_newer_file(blastfile, closuredir):
        mkdir(closuredir, overwrite=True)
        closure = True

    if closure:
        idsfile = op.join(closuredir, prefix + ".ids")
        cmd = "cut -f2 {0} | sort -u".format(blastfile)
        sh(cmd, outfile=idsfile)

        idsfastafile = op.join(closuredir, prefix + ".ids.fasta")
        cmd = "faSomeRecords {0} {1} {2}".format(poolfasta, idsfile, idsfastafile)
        sh(cmd)

        mergedfastafile = op.join(closuredir, prefix + ".merged.fasta")
        cmd = "cat {0} {1}".format(ctgfasta, idsfastafile)
        sh(cmd, outfile=mergedfastafile)

        afgfile = op.join(closuredir, prefix + ".afg")
        cmd = "toAmos -s {0} -o {1}".format(mergedfastafile, afgfile)
        sh(cmd)

        cwd = os.getcwd()
        os.chdir(closuredir)
        cmd = "minimus2 {0} -D REFCOUNT=0".format(prefix)
        sh(cmd)
        os.chdir(cwd)

    # Analyze output, make sure that:
    # + Get the singletons of the original set back
    # + Drop any contig that is comprised entirely of pulled set
    originalIDs = set(Fasta(ctgfasta).iterkeys())
    minimuscontig = op.join(closuredir, prefix + ".contig")
    c = ContigFile(minimuscontig)
    excludecontigs = set()
    for rec in c.iter_records():
        reads = set(x.id for x in rec.reads)
        if reads.isdisjoint(originalIDs):
            excludecontigs.add(rec.id)

    logging.debug("Exclude contigs: {0}".\
            format(", ".join(sorted(excludecontigs))))

    finalfasta = prefix + ".improved.fasta_"
    fw = open(finalfasta, "w")
    minimusfasta = op.join(closuredir, prefix + ".fasta")
    f = Fasta(minimusfasta)
    for id, rec in f.iteritems_ordered():
        if id in excludecontigs:
            continue
        SeqIO.write([rec], fw, "fasta")

    singletonfile = op.join(closuredir, prefix + ".singletons")
    singletons = set(x.strip() for x in open(singletonfile))
    leftovers = singletons & originalIDs

    logging.debug("Pull leftover singletons: {0}".\
            format(", ".join(sorted(leftovers))))

    f = Fasta(ctgfasta)
    for id, rec in f.iteritems_ordered():
        if id not in leftovers:
            continue
        SeqIO.write([rec], fw, "fasta")

    fw.close()

    fastafile = finalfasta
    finalfasta = fastafile.rstrip("_")
    format([fastafile, finalfasta, "--sequential", "--pad0=3",
        "--prefix={0}_".format(rid)])

    logging.debug("Improved FASTA written to `{0}`.".format(finalfasta))

    n50([ctgfasta])
    n50([finalfasta])

    os.remove(fastafile)
    os.remove(blastfile)
    os.remove("error.log")


if __name__ == '__main__':
    main()
