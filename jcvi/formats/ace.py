#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Parse ace format files.
"""
from __future__ import print_function

import os
import sys
import logging

from Bio.Sequencing import Ace
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcvi.apps.base import OptionParser, ActionDispatcher
from jcvi.formats.base import must_open


def main():

    actions = (
        ("extract", "extract contigs from ace file"),
        (
            "report",
            "generate report of read location, consensus location or quality segment per contig",
        ),
    )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def extract(args):
    """
    %prog extract [--options] ace_file

    Extract contigs from ace file and if necessary reformat header with
    a pipe(|) separated list of constituent reads.
    """
    p = OptionParser(extract.__doc__)
    p.add_option(
        "--format",
        default=False,
        action="store_true",
        help="enable flag to reformat header into a symbol separated list of constituent reads "
        + "[default: %default]",
    )
    p.add_option(
        "--singlets",
        default=False,
        action="store_true",
        help="ask the program to look in the singlets file (should be in the same folder) for "
        + "unused reads and put them in the resultant fasta file",
    )
    p.set_sep(sep="|", help="Separator used to list the reads in the FASTA header")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (acefile,) = args
    ace = Ace.read(must_open(acefile))
    logging.debug("Loaded ace file {0}".format(acefile))

    fastafile = acefile.rsplit(".", 1)[0] + ".fasta"
    fw = open(fastafile, "w")
    for c in ace.contigs:
        id = c.name
        if opts.format:
            id = opts.sep.join([read.name for read in c.af])

        seqrec = SeqRecord(Seq(c.sequence), id=id, description="")
        SeqIO.write([seqrec], fw, "fasta")

    if opts.singlets:
        singletsfile = acefile.rsplit(".", 1)[0] + ".singlets"
        if os.path.getsize(singletsfile) > 0:
            fp = SeqIO.parse(must_open(singletsfile), "fasta")
            for rec in fp:
                SeqIO.write(rec, fw, "fasta")

    fw.close()
    logging.debug("Wrote contigs to fasta file {0}".format(fastafile))


def report(args):
    """
    %prog report [--options] ace_file > report

    Prepare a report of read location, consensus location or quality segment per contig
    """
    from jcvi.utils.table import tabulate

    p = OptionParser(report.__doc__)

    types = {
        "read": ["padded_start", "padded_end", "orient"],
        "consensus": ["padded_consensus_start", "padded_consensus_end"],
        "quality": [
            "qual_clipping_start",
            "qual_clipping_end",
            "align_clipping_start",
            "align_clipping_end",
        ],
    }
    valid_types = tuple(types.keys())
    p.add_option(
        "--type", default="read", choices=valid_types, help="choose report type",
    )

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (acefile,) = args
    ace = Ace.read(must_open(acefile))
    logging.debug("Loaded ace file {0}".format(acefile))

    for c in ace.contigs:
        print(c.name)
        table = dict()
        if opts.type == "read":
            ps, pe = [], []
            ps = [read.padded_start for read in c.af]
            for i in range(1, len(ps)):
                pe.append(ps[i] - ps[i - 1])
            pe.append(c.nbases)
            map = dict(zip(ps, pe))
            for i, read in enumerate(c.af):
                values = [
                    str(x)
                    for x in (read.padded_start, map[read.padded_start], read.coru)
                ]
                for i, label in enumerate(types[opts.type]):
                    table[(str(read.name), label)] = values[i]
        elif opts.type == "consensus":
            for read in c.bs:
                values = [str(x) for x in (read.padded_start, read.padded_end)]
                for i, label in enumerate(types[opts.type]):
                    table[(str(read.name), label)] = values[i]
        elif opts.type == "quality":
            for read in c.reads:
                (r1, r2) = (read.rd, read.qa)
                values = [
                    str(x)
                    for x in (
                        r2.qual_clipping_start,
                        r2.qual_clipping_end,
                        r2.align_clipping_start,
                        r2.align_clipping_end,
                    )
                ]
                for i, label in enumerate(types[opts.type]):
                    table[(str(r1.name), label)] = values[i]
        print(tabulate(table), "\n")
    # logging.debug("Table with {0} rows written to `{1}`.".format(nrows, csvfile))


if __name__ == "__main__":
    main()
