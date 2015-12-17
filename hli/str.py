#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Identify repeat numbers in STR repeats.
"""

import os.path as op
import sys
import logging
import pysam

from collections import defaultdict
from jcvi.graphics.histogram import stem_leaf_plot
from jcvi.utils.cbook import percentage
from jcvi.apps.grid import MakeManager
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir


def main():

    actions = (
        ('pe', 'infer paired-end reads spanning a certain region'),
        ('lobstrindex', 'make lobSTR index'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def lobstrindex(args):
    """
    %prog lobstrindex hg38.trf.bed hg38.upper.fa

    Make lobSTR index. Make sure the FASTA contain only upper case (so use
    fasta.format --upper to convert from UCSC fasta).
    """
    p = OptionParser(lobstrindex.__doc__)
    p.set_home("lobstr")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    trfbed, fastafile = args
    lhome = opts.lobstr_home
    pf = fastafile.split(".")[0]
    mkdir(pf)

    newbedfile = trfbed + ".new"
    newbed = open(newbedfile, "w")
    fp = open(trfbed)
    retained = total = 0
    for row in fp:
        atoms = row.split()
        start = int(atoms[1]) + 1
        period = int(atoms[4])
        copynum = int(float(atoms[5]))
        length = period * copynum
        total += 1
        if period > 6 or length > 150:  # lobSTR cannot type these loci
            continue
        retained += 1
        end = start + length - 1
        print >> newbed, "\t".join(str(x) for x in (atoms[0], start, end,
                                period, copynum, "\t".join(atoms[6:])))
    newbed.close()
    logging.debug("Retained: {0}".format(percentage(retained, total)))

    mm = MakeManager()
    cmd = "python {0}/scripts/lobstr_index.py".format(lhome)
    cmd += " --str {0} --ref {1} --out {2}".format(newbedfile, fastafile, pf)
    mm.add((newbedfile, fastafile), op.join(pf, "lobSTR_ref.fasta.rsa"), cmd)

    tabfile = "{0}/{0}.tab".format(pf)
    cmd = "python {0}/scripts/GetSTRInfo.py".format(lhome)
    cmd += " {0} {1} > {2}".format(newbedfile, fastafile, tabfile)
    mm.add((newbedfile, fastafile), tabfile, cmd)
    mm.write()


def pe(args):
    """
    %prog pe bam chr start end

    Infer distance paired-end reads spanning a certain region.
    """
    p = OptionParser(pe.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    bam, chr, start, end = args
    start, end = int(start), int(end)
    target_len = end - start + 1
    pad = 1000  # How far do we look beyond the target
    pstart = start - pad
    pend = end + pad
    logging.debug("Target length={0}".format(target_len))

    samfile = pysam.AlignmentFile(bam, "rb")
    iter = samfile.fetch(chr, pstart, pend)
    cache = defaultdict(list)
    for x in iter:
        if not x.is_paired:
            continue
        cache[x.query_name].append(x)

    tlens = []
    for name, reads in cache.iteritems():
        if len(reads) < 2:
            continue
        a, b = reads[:2]
        # Get all pairs where read1 is on left flank and read2 is on right flank
        if a.reference_start >= start or b.reference_end <= end:
            continue
        for x in (a, b):
            print x.query_name, x.is_reverse, \
                    x.query_alignment_start, x.query_alignment_end, \
                    x.reference_start, x.reference_end, \
                    x.tlen
        print '=' * 60
        tlen = abs(x.tlen)
        tlens.append(tlen)

    stem_leaf_plot(tlens, 300, 600, 20)


if __name__ == '__main__':
    main()
