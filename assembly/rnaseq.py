#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import sys

from collections import defaultdict
from jcvi.apps.base import MOptionParser

from jcvi.formats.blast import Blast
from jcvi.formats.sizes import Sizes
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('qc', 'evaluate completeness, contiguity, and chimer of RNAseq'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def qc(args):
    """
    %prog rnaseq blastfile ref.fasta

    Evaluate de-novo RNA-seq assembly against a reference gene set (same or
    closely related organism). Ideally blatfile needs to be supermap'd.

    Following metric is used (Martin et al. 2010, Rnnotator paper):
    Accuracy: % of contigs share >=95% identity with ref genome (TODO)
    Completeness: % of ref genes covered by contigs to >=80% of their lengths
    Contiguity: % of ref genes covered by a *single* contig >=80% of lengths
    Chimer: % of contigs that contain two or more annotated genes >= 50bp
    """
    from jcvi.algorithms.supermap import supermap

    p = MOptionParser(qc.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    blastfile, reffasta = args
    sizes = Sizes(reffasta).mapping
    known_genes = len(sizes)

    querysupermap = blastfile + ".query.supermap"
    refsupermap = blastfile + ".ref.supermap"

    if not op.exists(querysupermap):
        supermap(blastfile, filter="query")
    if not op.exists(refsupermap):
        supermap(blastfile, filter="ref")

    blast = Blast(querysupermap)
    chimers = 0
    goodctg80 = set()
    goodctg50 = set()
    for ctg, hits in blast.iter_hits():
        bps = defaultdict(int)
        for x in hits:
            bps[x.subject] += abs(x.sstop - x.sstart) + 1

        valid_hits = bps.items()
        for vh, length in valid_hits:
            rsize = sizes[vh]
            ratio = length * 100. / rsize
            if ratio >= 80:
                goodctg80.add(ctg)
            if ratio >= 50:
                goodctg50.add(ctg)

        # Chimer
        if len(valid_hits) > 1:
            chimers += 1

    blast = Blast(refsupermap)
    goodref80 = set()
    goodref50 = set()
    bps = defaultdict(int)
    for x in blast:
        bps[x.subject] += abs(x.sstop - x.sstart) + 1

    for vh, length in bps.items():
        rsize = sizes[vh]
        ratio = length * 100. / rsize
        if ratio >= 80:
            goodref80.add(vh)
        if ratio >= 50:
            goodref50.add(vh)

    print >> sys.stderr, "Reference set: `{0}`,  # of transcripts {1}".\
            format(reffasta, known_genes)
    print >> sys.stderr, "A total of {0} contigs map to 80% of a reference"\
            " transcript".format(len(goodctg80))
    print >> sys.stderr, "A total of {0} contigs map to 50% of a reference"\
            " transcript".format(len(goodctg50))
    print >> sys.stderr, "A total of {0} reference transcripts ({1:.1f}%) have 80% covered" \
            .format(len(goodref80), len(goodref80) * 100. / known_genes)


if __name__ == '__main__':
    main()
