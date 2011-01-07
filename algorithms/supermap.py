#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] infile
This script combines a series of pairwise alignments, sort and filter the alignments

In order to handle dups, we have to run two monotonous chains in both genomes,
first chain using ref1, and a second chain using ref2 and we will have options
to keep either the union or the intersection of retain chained alignments from
both genomes, similar to the SUPERMAP algorithm. This operation is symmetrical.

Finally a filtered subset will be written to outfile.
"""

import os
import sys
import collections
import logging

from optparse import OptionParser

from bx.align import maf
from bx import interval_index_file

Weighted_interval = collections.namedtuple("Weighted_interval",
        "chromosome left right weight")

# Stolen from Brad Chapman's blog
# http://bcbio.wordpress.com/2009/07/26/sorting-genomic-alignments-using-python/
def build_index(in_file, index_file):
    indexes = interval_index_file.Indexes()
    with open(in_file) as in_handle:
        reader = maf.Reader(in_handle)
        while 1:
            pos = reader.file.tell()
            rec = reader.next()
            if rec is None:
                break
            for c in rec.components:
                indexes.add(c.src, c.forward_strand_start,
                        c.forward_strand_end, pos, max=c.src_size )

    with open(index_file, "w") as index_handle:
        indexes.write(index_handle)


def interval_chain(intervals, endpoints):
    """
    Take a list of weighted intervals, return a non-overlapping set with max weight
    We proceed by looking at the each end point (sorted by their relative positions)
    """
    endpoints.sort()
    left_index = {}   # stores the left end index for quick retrieval

    # dynamic programming, each entry [score, from_index, which_chain]
    scores = [[0, -1, -1]]
    for i, e in enumerate(endpoints):
        pos, j, weight = e  # endpoint bp position, interval id, weight
        cur_score = scores[-1][:]

        if weight > 0: # this is right end of j-th interval
            # update if chaining j-th interval gives a better score
            left_j = left_index[j]
            chain_score = scores[left_j][0] + weight
            if chain_score > cur_score[0]: 
                cur_score = [chain_score, left_j, j]
        else:  # left end
            left_index[j] = i+1
        scores.append(cur_score)

    chains = []
    score, last, chain_id = scores[-1] # start backtracking
    while last!=-1:
        if chain_id!=-1: 
            chains.append(chain_id)
        _, last, chain_id = scores[last]

    chains.reverse()
    print len(chains), "intervals retained with score %d" % score
    total_aligned_bases = 0
    for c in chains:
        #print intervals[c]
        total_aligned_bases += intervals[c].right - intervals[c].left + 1
    print "total aligned bases are %d" % total_aligned_bases

    return set([c/2 for c in chains]) # remember consecutive chains belong to the
                                      # same alignment


def main(options, args):
    
    in_file = args[0]
    base, ext = os.path.splitext(in_file)
    out_file = "%s-filtered%s" %(base, ext)
    index_file = in_file + ".index"
    if not os.path.exists(index_file):
        build_index(in_file, index_file)
    index = maf.Indexed(in_file, index_file)

    fp = file(in_file)
    reader = maf.Reader(fp)

    intervals = [] # give each interval a unique id
    endpoints = collections.defaultdict(list) # chromosome => list of endpoints 
    filtered_rec = set()
    j = 0
    rec_info = []
    while 1:
        pos = reader.file.tell()
        rec_info.append((j/2, pos))   # position of alignment j in file
        rec = reader.next()
        if rec is None:
            break
        for c in rec.components:
            chromosome, left, right, weight = c.src, c.forward_strand_start, \
                    c.forward_strand_end, rec.score

            interval = Weighted_interval(chromosome, left, right, weight)
            intervals.append(interval)
            endpoints[chromosome].append((left, j, -weight))  # left end
            endpoints[chromosome].append((right, j, weight))  # right end
            j += 1

    for chromosome in sorted(endpoints.keys()):
        v = endpoints[chromosome]
        print chromosome, ": start with %d intervals" % (len(v)/2) 
        filtered_rec |= interval_chain(intervals, v)

    print "filtered alignment size %d" % len(filtered_rec)

    fw = file(out_file, "w")
    writer = maf.Writer(fw)

    for j, pos in rec_info:
        if j in filtered_rec:
            rec = index.get_at_offset(pos)
            writer.write(rec)

    fp.close()
    fw.close()


if __name__ == '__main__':
    
    p = OptionParser(__doc__)

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(parser.print_help())

    main(options, args)
