#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import logging
import collections

import numpy as np
from optparse import OptionParser

from jcvi.formats.bed import Bed
from jcvi.formats.blast import BlastLine
from jcvi.utils.grouper import Grouper
from jcvi.apps.base import ActionDispatcher, debug
debug()


def _score(cluster):
    """ 
    score of the cluster, in this case, is the number of non-repetitive matches
    """
    x, y = zip(*cluster)
    return min(len(set(x)), len(set(y)))


def group_hits(blasts):
    # grouping the hits based on chromosome pair
    all_hits = collections.defaultdict(list)
    for b in blasts:
        all_hits[(b.qseqid, b.sseqid)].append((b.qi, b.si))

    return all_hits


def read_blast(blast_file, qorder, sorder, is_self=False):
    """
    read the blast and convert name into coordinates
    """
    fp = open(blast_file)
    filtered_blast = []
    seen = set()
    for row in fp:
        b = BlastLine(row)
        query, subject = b.query, b.subject
        if query not in qorder or subject not in sorder: continue

        key = query, subject
        if key in seen: continue
        seen.add(key)

        qi, q = qorder[query]
        si, s = sorder[subject]

        if is_self and qi > si:
            # remove redundant a<->b to one side when doing self-self BLAST
            query, subject = subject, query
            qi, si = si, qi
            q, s = s, q

        b.qseqid, b.sseqid = q.seqid, s.seqid
        b.qi, b.si = qi, si

        filtered_blast.append(b)

    return filtered_blast


def read_anchors(anchor_file, qorder, sorder):
    """
    anchors file are just (geneA, geneB) pairs (with possible deflines)
    """
    all_anchors = collections.defaultdict(list)
    fp = open(anchor_file)
    for row in fp:
        if row[0]=='#': continue
        a, b = row.split()
        if a not in qorder or b not in sorder: continue
        qi, q = qorder[a]
        si, s = sorder[b]
        all_anchors[(q.seqid, s.seqid)].append((qi, si))

    return all_anchors


def synteny_scan(points, xdist, ydist, N):
    """
    This is the core single linkage algorithm
    this behaves in O(n) complexity: we iterate through the pairs, for each pair
    we look back on the adjacent pairs to find links
    """
    clusters = Grouper()
    n = len(points)
    points.sort()
    for i in xrange(n):
        for j in xrange(i-1, -1, -1):
            # x-axis distance
            del_x = points[i][0]-points[j][0]
            if del_x > xdist: break
            # y-axis distance
            del_y = points[i][1]-points[j][1]
            if abs(del_y) > ydist: continue
            # otherwise join
            clusters.join(points[i], points[j])

    # select clusters that are at least >=N
    clusters = [sorted(cluster) for cluster in list(clusters) \
            if _score(cluster)>=N]

    return clusters


def batch_scan(points, xdist=20, ydist=20, N=6):
    """
    runs synteny_scan() per chromosome pair
    """
    chr_pair_points = group_hits(points)

    clusters = []
    for chr_pair in sorted(chr_pair_points.keys()):
        points = chr_pair_points[chr_pair]
        logging.debug("%s: %d" % (chr_pair, len(points)))
        clusters.extend(synteny_scan(points, xdist, ydist, N))

    return clusters


def synteny_liftover(points, anchors, dist):
    """
    This is to get the nearest anchors for all the points (useful for the
    `liftover` operation below).
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        raise ImportError("You must install python package `scipy` " + \
                "(http://www.scipy.org)")

    points = np.array(points)
    ppoints = points[:, :2] if points.shape[1] > 2 else points
    tree = cKDTree(anchors, leafsize=16)
    #print tree.data
    dists, idxs = tree.query(ppoints, p=1, distance_upper_bound=dist)
    #print [(d, idx) for (d, idx) in zip(dists, idxs) if idx!=tree.n]

    for point, dist, idx in zip(points, dists, idxs):
        # nearest is out of range
        if idx==tree.n: continue
        yield point


def add_options(p, args):
    """
    scan and liftover has similar interfaces, so share common options
    returns opts, files
    """
    p.add_option("--qbed", dest="qbed", help="path to qbed (required)")
    p.add_option("--sbed", dest="sbed", help="path to sbed (required)")

    p.add_option("--dist", dest="dist",
            default=10, type="int", 
            help="the extent of flanking regions to search [default: %default]")

    opts, args = p.parse_args(args)

    if not (len(args) == 2 and opts.qbed and opts.sbed):
        sys.exit(p.print_help())

    blast_file, anchor_file = args 

    qbed_file, sbed_file = opts.qbed, opts.sbed
    # is this a self-self blast?
    is_self = (qbed_file == sbed_file)
    if is_self:
        logging.debug("Looks like self-self BLAST")

    return blast_file, anchor_file, qbed_file, sbed_file, \
            opts.dist, is_self, opts


def main():

    actions = (
        ('scan', 'get anchor list using single-linkage algorithm'),
        ('liftover', 'given anchor list, pull adjancent pairs from blast file')
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def scan(args):
    """
    %prog scan blastfile anchor_file [options]

    pull out syntenic anchors from blastfile based on single-linkage algorithm
    """
    p = OptionParser(scan.__doc__)
    p.add_option("-N", dest="N", type="int", default=5,
            help="minimum number of anchors in a cluster [default: %default]")

    blast_file, anchor_file, qbed_file, sbed_file, dist, is_self, opts = \
            add_options(p, args)

    qbed = Bed(qbed_file)
    sbed = Bed(sbed_file)
    qorder = qbed.order
    sorder = sbed.order

    filtered_blast = read_blast(blast_file, qorder, sorder, is_self=is_self)

    fw = open(anchor_file, "w")
    clusters = batch_scan(filtered_blast, xdist=dist, ydist=dist, N=opts.N)
    for cluster in clusters:
        print >>fw, "###"
        for qi, si in cluster:
            query, subject = qbed[qi].accn, sbed[si].accn
            print >>fw, "\t".join((query, subject))


def liftover(args):
    """
    %prog liftover blastfile anchorfile [options] 

    Typical use for this program is given a list of anchors (for example, syntennic
    genes), choose from the blastfile the pairs that are close to the anchors.

    Anchorfile has the following format, each row defines a pair.

        geneA geneB
        geneC geneD
    """
    p = OptionParser(liftover.__doc__)

    blast_file, anchor_file, qbed_file, sbed_file, dist, is_self, opts = \
            add_options(p, args)

    qbed = Bed(qbed_file)
    sbed = Bed(sbed_file)
    qorder = qbed.order
    sorder = sbed.order

    filtered_blast = read_blast(blast_file, qorder, sorder, is_self=is_self)
    all_hits = group_hits(filtered_blast)
    all_anchors = read_anchors(anchor_file, qorder, sorder)

    # select hits that are close to the anchor list
    j = 0
    fw = sys.stdout
    for chr_pair in sorted(all_anchors.keys()):
        hits = np.array(all_hits[chr_pair])
        anchors = np.array(all_anchors[chr_pair])

        logging.debug("%s: %d" % (chr_pair, len(anchors)))
        if not len(hits): continue

        for point in synteny_liftover(hits, anchors, dist):
            qi, si = point[:2]
            query, subject = qbed[qi].accn, sbed[si].accn
            print >>fw, "\t".join((query, subject, "lifted"))
            j+=1
    
    logging.debug("%d new pairs found" % j)


if __name__ == '__main__':
    main()
    
