#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run EMBOSS programs.
"""

import sys

from jcvi.apps.base import OptionParser, ActionDispatcher


class NeedleHeader (object):

    def __init__(self, filename):
        fp = open(filename)
        for row in fp:
            if row[0] != '#':
                continue
            # Identity:      89/89 (100.0%)
            if row.startswith('# Identity'):
                self.identity = row.split(":")[-1].strip()
            if row.startswith('# Score'):
                self.score = row.split(":")[-1].strip()


def main():

    actions = (
        ('needle', 'take protein pairs and needle them'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def needle(args):
    """
    %prog needle pairs a.pep.fasta b.pep.fasta

    Take protein pairs and needle them.
    """
    from Bio.Emboss.Applications import NeedleCommandline

    from jcvi.formats.fasta import Fasta, SeqIO
    from jcvi.formats.base import FileShredder

    p = OptionParser(needle.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    pairsfile, apep, bpep = args
    afasta = Fasta(apep)
    bfasta = Fasta(bpep)
    fp = open(pairsfile)
    for row in fp:
        fa = open(pairsfile + "_a.fasta", "w")
        fb = open(pairsfile + "_b.fasta", "w")
        a, b = row.split()
        a = afasta[a]
        b = bfasta[b]
        SeqIO.write([a], fa, "fasta")
        SeqIO.write([b], fb, "fasta")
        fa.close()
        fb.close()
        needlefile = pairsfile + "_ab.needle"
        needle_cline = NeedleCommandline(asequence=fa.name,
                            bsequence=fb.name,
                            gapopen=10, gapextend=0.5,
                            outfile=needlefile)
        stdout, stderr = needle_cline()
        print >> sys.stderr, stdout + stderr
        #align = AlignIO.read(needlefile, "emboss")
        nh = NeedleHeader(needlefile)
        print "\t".join((a.id, b.id, nh.identity, nh.score))
        FileShredder([fa.name, fb.name, needlefile])


if __name__ == '__main__':
    main()
