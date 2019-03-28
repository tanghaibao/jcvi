#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run EMBOSS programs.
"""
from __future__ import print_function

import sys
import multiprocessing as mp

from jcvi.apps.base import OptionParser, ActionDispatcher
from jcvi.formats.base import FileShredder, must_open

class NeedleHeader (object):

    def __init__(self, filename):
        fp = must_open(filename)
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


def _needle(fa, fb, needlefile, a, b, results):
    """
    Run single needle job
    """
    from Bio.Emboss.Applications import NeedleCommandline

    needle_cline = NeedleCommandline(asequence=fa, bsequence=fb,
                        gapopen=10, gapextend=0.5, outfile=needlefile)
    stdout, stderr = needle_cline()
    nh = NeedleHeader(needlefile)
    FileShredder([fa, fb, needlefile], verbose=False)
    r = ["\t".join((a, b, nh.identity, nh.score))]

    results.extend(r)


def needle(args):
    """
    %prog needle nw.pairs a.pep.fasta b.pep.fasta

    Take protein pairs and needle them
    Automatically writes output file `nw.scores`
    """
    from jcvi.formats.fasta import Fasta, SeqIO

    p = OptionParser(needle.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    manager = mp.Manager()
    results = manager.list()
    needle_pool = mp.Pool(processes=mp.cpu_count())

    pairsfile, apep, bpep = args
    afasta, bfasta = Fasta(apep), Fasta(bpep)
    fp = must_open(pairsfile)
    for i, row in enumerate(fp):
        a, b = row.split()
        a, b = afasta[a], bfasta[b]
        fa, fb = must_open("{0}_{1}_a.fasta".format(pairsfile, i), "w"), \
            must_open("{0}_{1}_b.fasta".format(pairsfile, i), "w")
        SeqIO.write([a], fa, "fasta")
        SeqIO.write([b], fb, "fasta")
        fa.close()
        fb.close()

        needlefile = "{0}_{1}_ab.needle".format(pairsfile, i)
        needle_pool.apply_async(_needle, \
            (fa.name, fb.name, needlefile, a.id, b.id, results))

    needle_pool.close()
    needle_pool.join()

    fp.close()

    scoresfile = "{0}.scores".format(pairsfile.rsplit(".")[0])
    fw = must_open(scoresfile, "w")
    for result in results:
        print(result, file=fw)
    fw.close()


if __name__ == '__main__':
    main()
