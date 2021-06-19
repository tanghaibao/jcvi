#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This implements the method described in Tang et al. 2010 PNAS paper,
<http://www.pnas.org/content/107/1/472>

Angiosperm genome comparisons reveal early polyploidy in the monocot lineage

The main pipeline assumes starting with defined synteny blocks in .anchors
format (use compara.synteny.scan()), then segment the chromosomes and cluster
segments according to the matching patterns. Finally the putative ancestral
regions (PAR) are identified and visualized.
"""
import os.path as op
import sys
import logging
from math import log

import numpy as np
from more_itertools import pairwise

from jcvi.compara.synteny import AnchorFile, check_beds
from jcvi.formats.bed import Bed
from jcvi.formats.blast import BlastLine
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh


def main():

    actions = (
        ("cluster", "cluster the segments"),
        ("pad", "test and reconstruct candidate PADs"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def make_arrays(blastfile, qpadbed, spadbed, qpadnames, spadnames):
    """
    This function makes three matrices: observed, expected and logmp. The logmp
    contains the statistical significance for each comparison.
    """
    m, n = len(qpadnames), len(spadnames)
    qpadorder, spadorder = qpadbed.order, spadbed.order
    qpadid = dict((a, i) for i, a in enumerate(qpadnames))
    spadid = dict((a, i) for i, a in enumerate(spadnames))
    qpadlen = dict((a, len(b)) for a, b in qpadbed.sub_beds())
    spadlen = dict((a, len(b)) for a, b in spadbed.sub_beds())

    qsize, ssize = len(qpadbed), len(spadbed)

    assert sum(qpadlen.values()) == qsize
    assert sum(spadlen.values()) == ssize

    # Populate arrays of observed counts and expected counts
    logging.debug("Initialize array of size ({0} x {1})".format(m, n))
    observed = np.zeros((m, n))
    fp = open(blastfile)
    all_dots = 0
    for row in fp:
        b = BlastLine(row)
        qi, q = qpadorder[b.query]
        si, s = spadorder[b.subject]
        qseqid, sseqid = q.seqid, s.seqid
        qsi, ssi = qpadid[qseqid], spadid[sseqid]
        observed[qsi, ssi] += 1
        all_dots += 1

    assert int(round(observed.sum())) == all_dots

    logging.debug("Total area: {0} x {1}".format(qsize, ssize))
    S = qsize * ssize
    expected = np.zeros((m, n))
    qsum = 0
    for i, a in enumerate(qpadnames):
        alen = qpadlen[a]
        qsum += alen
        for j, b in enumerate(spadnames):
            blen = spadlen[b]
            expected[i, j] = all_dots * alen * blen * 1.0 / S

    assert int(round(expected.sum())) == all_dots

    # Calculate the statistical significance for each cell
    from scipy.stats.distributions import poisson

    logmp = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            obs, exp = observed[i, j], expected[i, j]
            pois = max(poisson.pmf(obs, exp), 1e-250)  # Underflow
            logmp[i, j] = max(-log(pois), 0)

    return logmp


def pad(args):
    """
    %prog pad blastfile cdtfile --qbed q.pad.bed --sbed s.pad.bed

    Test and reconstruct candidate PADs.
    """
    from jcvi.formats.cdt import CDT

    p = OptionParser(pad.__doc__)
    p.set_beds()
    p.add_option(
        "--cutoff",
        default=0.3,
        type="float",
        help="The clustering cutoff to call similar",
    )

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    cutoff = opts.cutoff
    blastfile, cdtfile = args
    qbed, sbed, qorder, sorder, is_self = check_beds(blastfile, p, opts)

    cdt = CDT(cdtfile)
    qparts = list(cdt.iter_partitions(cutoff=cutoff))
    sparts = list(cdt.iter_partitions(cutoff=cutoff, gtr=False))

    qid, sid = {}, {}
    for i, part in enumerate(qparts):
        qid.update(dict((x, i) for x in part))
    for i, part in enumerate(sparts):
        sid.update(dict((x, i) for x in part))

    # Without writing files, conversion from PAD to merged PAD is done in memory
    for q in qbed:
        q.seqid = qid[q.seqid]
    for s in sbed:
        s.seqid = sid[s.seqid]

    qnames = range(len(qparts))
    snames = range(len(sparts))

    logmp = make_arrays(blastfile, qbed, sbed, qnames, snames)
    m, n = logmp.shape
    pvalue_cutoff = 1e-30
    cutoff = -log(pvalue_cutoff)

    significant = []
    for i in range(m):
        for j in range(n):
            score = logmp[i, j]
            if score < cutoff:
                continue
            significant.append((qparts[i], sparts[j], score))

    for a, b, score in significant:
        print("|".join(a), "|".join(b), score)

    logging.debug(
        "Collected {0} PAR comparisons significant at (P < {1}).".format(
            len(significant), pvalue_cutoff
        )
    )

    return significant


def get_segments(ranges, extra, minsegment=40):
    """
    Given a list of Range, perform chaining on the ranges and select a highest
    scoring subset and cut based on their boundaries. Let's say the projection
    of the synteny blocks onto one axis look like the following.

    1=====10......20====30....35====~~

    Then the segmentation will yield a block [1, 20), [20, 35), using an
    arbitrary right extension rule. Extra are additional end breaks for
    chromosomes.
    """
    from jcvi.utils.range import range_chain, LEFT, RIGHT

    NUL = 2
    selected, score = range_chain(ranges)

    endpoints = [(x.start, NUL) for x in selected]
    endpoints += [(x[0], LEFT) for x in extra]
    endpoints += [(x[1], RIGHT) for x in extra]
    endpoints.sort()

    current_left = 0
    for a, ai in endpoints:

        if ai == LEFT:
            current_left = a
        if ai == RIGHT:
            yield current_left, a
        elif ai == NUL:
            if a - current_left < minsegment:
                continue
            yield current_left, a - 1
            current_left = a


def write_PAD_bed(bedfile, prefix, pads, bed):

    fw = open(bedfile, "w")
    padnames = ["{0}:{1:05d}-{2:05d}".format(prefix, a, b) for a, b in pads]
    for a, b in pairwise(padnames):
        assert a != b, a

    j = 0
    # Assign all genes to new partitions
    for i, x in enumerate(bed):
        if i > b:
            j += 1
        print("\t".join((padnames[j], str(i), str(i + 1), x.accn)), file=fw)

    fw.close()

    npads = len(pads)
    logging.debug("{0} partition written in `{1}`.".format(npads, bedfile))
    return npads, padnames


def cluster(args):
    """
    %prog cluster blastfile anchorfile --qbed qbedfile --sbed sbedfile

    Cluster the segments and form PAD. This is the method described in Tang et
    al. (2010) PNAS paper. The anchorfile defines a list of synteny blocks,
    based on which the genome on one or both axis can be chopped up into pieces
    and clustered.
    """
    from jcvi.utils.range import Range

    p = OptionParser(cluster.__doc__)
    p.set_beds()
    p.add_option(
        "--minsize", default=10, type="int", help="Only segment using blocks >= size"
    )
    p.add_option(
        "--path", default="~/scratch/bin", help="Path to the CLUSTER 3.0 binary"
    )

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blastfile, anchorfile = args
    qbed, sbed, qorder, sorder, is_self = check_beds(blastfile, p, opts)

    minsize = opts.minsize
    ac = AnchorFile(anchorfile)
    qranges, sranges = [], []
    qextra = [x[1:] for x in qbed.get_breaks()]
    sextra = [x[1:] for x in sbed.get_breaks()]

    id = 0
    for block in ac.iter_blocks(minsize=minsize):
        q, s = list(zip(*block))[:2]
        q = [qorder[x][0] for x in q]
        s = [sorder[x][0] for x in s]
        minq, maxq = min(q), max(q)
        mins, maxs = min(s), max(s)
        id += 1

        qr = Range("0", minq, maxq, maxq - minq, id)
        sr = Range("0", mins, maxs, maxs - mins, id)
        qranges.append(qr)
        sranges.append(sr)

    qpads = list(get_segments(qranges, qextra))
    spads = list(get_segments(sranges, sextra))

    suffix = ".pad.bed"
    qpf = opts.qbed.split(".")[0]
    spf = opts.sbed.split(".")[0]
    qpadfile = qpf + suffix
    spadfile = spf + suffix
    qnpads, qpadnames = write_PAD_bed(qpadfile, qpf, qpads, qbed)
    snpads, spadnames = write_PAD_bed(spadfile, spf, spads, sbed)

    qpadbed, spadbed = Bed(qpadfile), Bed(spadfile)

    logmp = make_arrays(blastfile, qpadbed, spadbed, qpadnames, spadnames)
    m, n = logmp.shape

    matrixfile = ".".join((qpf, spf, "logmp.txt"))
    fw = open(matrixfile, "w")
    header = ["o"] + spadnames
    print("\t".join(header), file=fw)
    for i in range(m):
        row = [qpadnames[i]] + ["{0:.1f}".format(x) for x in logmp[i, :]]
        print("\t".join(row), file=fw)

    fw.close()

    # Run CLUSTER 3.0 (Pearson correlation, average linkage)
    cmd = op.join(opts.path, "cluster")
    cmd += " -g 2 -e 2 -m a -f {0}".format(matrixfile)
    pf = matrixfile.rsplit(".", 1)[0]
    cdtfile = pf + ".cdt"
    if need_update(matrixfile, cdtfile):
        sh(cmd)


if __name__ == "__main__":
    main()
