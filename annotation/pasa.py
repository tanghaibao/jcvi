#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Utilities when processing PASA results.
"""

import sys
import logging

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('longest', 'label longest transcript per gene as full-length'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def longest(args):
    """
    %prog longest pasa.fasta output.subclusters.out

    Find the longest PASA assembly and label it as full-length. Also removes
    transcripts shorter than half the length of the longest, or shorter than
    200bp. The assemblies for the same locus is found in
    `output.subclusters.out`. In particular the lines that look like:

    sub-cluster: asmbl_25 asmbl_26 asmbl_27
    """
    from jcvi.formats.fasta import Fasta, SeqIO
    from jcvi.formats.sizes import Sizes

    p = OptionParser(longest.__doc__)
    p.add_option("--prefix", default="pasa",
                 help="Replace asmbl_ with prefix [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, subclusters = args
    prefix = fastafile.rsplit(".", 1)[0]

    idsfile = prefix + ".fl.ids"
    fw = open(idsfile, "w")
    sizes = Sizes(fastafile).mapping

    name_convert = lambda x: x.replace("asmbl", opts.prefix)

    keep = set()  # List of IDs to write
    fp = open(subclusters)
    nrecs = 0
    for row in fp:
        if not row.startswith("sub-cluster:"):
            continue
        asmbls = row.split()[1:]
        longest_asmbl = max(asmbls, key=lambda x: sizes[x])
        longest_size = sizes[longest_asmbl]
        print >> fw, name_convert(longest_asmbl)
        nrecs += 1
        cutoff = max(longest_size / 2, 200)
        keep.update(set(x for x in asmbls if sizes[x] >= cutoff))

    fw.close()
    logging.debug("{0} fl-cDNA records written to `{1}`.".format(nrecs, idsfile))

    f = Fasta(fastafile, lazy=True)
    newfastafile = prefix + ".clean.fasta"
    fw = open(newfastafile, "w")
    nrecs = 0
    for name, rec in f.iteritems_ordered():
        if name not in keep:
            continue

        rec.id = name_convert(name)
        rec.description = ""
        SeqIO.write([rec], fw, "fasta")
        nrecs += 1

    fw.close()
    logging.debug("{0} valid records written to `{1}`.".format(nrecs, newfastafile))


if __name__ == '__main__':
    main()
