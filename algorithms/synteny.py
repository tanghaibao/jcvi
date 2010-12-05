#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import logging
import numpy as np
import collections
from optparse import OptionParser

from jcvi.formats.bed import Bed
from jcvi.formats.blast import Blast
from jcvi.utils.grouper import Grouper
from jcvi.apps.base import ActionDispatcher


def score(cluster):
    """ 
    score of the cluster, in this case, is the number of non-repetitive matches
    """
    x, y = zip(*cluster)
    return min(len(set(x)), len(set(y)))


def single_linkage(points, xdist, ydist, N):
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
    clusters = [cluster for cluster in list(clusters) if score(cluster)>=N]
    return clusters


def main():

    actions = (
        ('scan', 'get anchor list using single-linkage algorithm'),
        ('liftover', 'given anchor list, pull adjancent pairs from blast file')
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def scan(args):
    """
    %prog scan blastfile anchorfile [options]

    pull out syntenic anchors from blastfile based on single-linkage algorithm
    """
    pass


def liftover(args):
    """
    %prog liftover blastfile anchorfile [options] 

    typical use for this program is given an anchor list (for example, syntennic
    genes), choose from the second blast_file for which the pairs are close to the anchors.

    Anchor list should have the following format, each row defines a pair:

        geneA geneB
        geneC geneD

    Use KD-tree for processing distance query.
    """
    p = OptionParser(liftover.__doc__)

    p.add_option("--qbed", dest="qbed", help="path to qbed")
    p.add_option("--sbed", dest="sbed", help="path to sbed")

    p.add_option("--dist", dest="dist",
            default=10, type="int", 
            help="the extent of flanking regions to search [default: %default]") 

    opts, files = p.parse_args(args)

    if not (len(files) == 2 and opts.qbed and opts.sbed):
        sys.exit(p.print_help())

    blast_file, anchor_file = files

    logging.basicConfig(level=logging.DEBUG)
    try:
        from scipy.spatial import cKDTree
    except ImportError, e:
        logging.error(e)
        logging.error("You must install python package `scipy` " + \
                "(http://www.scipy.org)")
        sys.exit(1)

    qbed_file, sbed_file = opts.qbed, opts.sbed
    # is this a self-self blast?
    is_self = (qbed_file == sbed_file)
    if is_self:
        logging.debug("... looks like a self-self BLAST to me")

    logging.debug("read annotation files %s and %s" % (qbed_file, sbed_file))
    qbed = Bed(qbed_file)
    sbed = Bed(sbed_file)
    qorder = qbed.get_order()
    sorder = sbed.get_order()

    blast = Blast(blast_file)
    filtered_blasts = []
    seen = set()
    for b in blast:
        #query, subject = _(b.query), _(b.subject)
        if query not in qorder or subject not in sorder: continue
        qi, q = qorder[query]
        si, s = sorder[subject]

        if is_self and qi > si:
            # remove redundant a<->b to one side when doing self-self BLAST
            query, subject = subject, query
            qi, si = si, qi
            q, s = s, q

        key = query, subject
        if key in seen: continue
        seen.add(key)
        b.query, b.subject = key

        b.qi, b.si = qi, si
        b.qseqid, b.sseqid = q.seqid, s.seqid
        
        filtered_blasts.append(b)

    all_anchors = collections.defaultdict(list)
    fp = file(anchor_file)
    for row in fp:
        if row[0]=='#': continue
        a, b = row.split()
        if a not in qorder or b not in sorder: continue
        qi, q = qorder[a]
        si, s = sorder[b]
        all_anchors[(q.seqid, s.seqid)].append((qi, si))

    # grouping the hits based on chromosome pair for sending in find_nearby
    all_hits = collections.defaultdict(list)
    for b in filtered_blasts:
        all_hits[(b.qseqid, b.sseqid)].append((b.qi, b.si))

    # select hits that are close to the anchor list
    j = 0
    fw = sys.stdout
    for chr_pair in sorted(all_hits.keys()):
        hits = np.array(all_hits[chr_pair])
        anchors = np.array(all_anchors[chr_pair])

        logging.debug("%s: %d" % (chr_pair, len(anchors)))
        if not anchors: continue
        tree = cKDTree(anchors, leafsize=16)
        #print tree.data
        dists, idxs = tree.query(hits, p=1, distance_upper_bound=opts.dist)
        #print [(d, idx) for (d, idx) in zip(dists, idxs) if idx!=tree.n]

        for i, (dd, idx) in enumerate(zip(dists, idxs)):
            if dd==0: continue # same anchors
            if idx!=tree.n:
                qi, si = hits[i]
                query, subject = qbed[qi]["accn"], sbed[si]["accn"]
                print >>fw, "\t".join((query, subject, "lifted"))
                j+=1
    
    logging.debug("%d new pairs found" % j)



if __name__ == '__main__':
    main()
    
