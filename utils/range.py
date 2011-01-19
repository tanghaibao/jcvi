#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This script implements algorithm for finding intersecting rectangles, 
both on the 2D dotplot and 1D-projection

`range_chain` implements the exon-chain algorithm
"""

from collections import namedtuple

from jcvi.utils.grouper import Grouper

LEFT, RIGHT = 0, 1 

Range = namedtuple("Range", "seqid start end score id")


def range_overlap(a, b):
    """
    >>> range_overlap(("1", 30, 45), ("1", 45, 55))
    True
    >>> range_overlap(("1", 30, 45), ("1", 57, 68))
    False
    >>> range_overlap(("1", 30, 45), ("2", 42, 55))
    False
    """
    a_chr, a_min, a_max = a
    b_chr, b_min, b_max = b
    # must be on the same chromosome
    if a_chr!=b_chr: return False 
    return (a_min <= b_max) and (b_min <= a_max)


def range_minmax(ranges):
    """
    Returns the span of a collection of ranges where start is the smallest of
    all starts, and end is the largest of all ends
    >>> ranges = [(30, 45), (40, 50), (10, 100)]
    >>> range_minmax(ranges)
    (10, 100)
    """
    rmin = min(ranges)[0]
    rmax = max(ranges, key=lambda x: x[1])[1]
    return rmin, rmax


def range_union(ranges):
    """
    Returns total size of ranges, expect range as (chr, left, right)
    >>> ranges = [("1", 30, 45), ("1", 40, 50), ("1", 10, 50)]
    >>> range_union(ranges)
    41
    >>> ranges = [("1", 30, 45), ("2", 40, 50)]
    >>> range_union(ranges)
    27
    >>> ranges = [("1", 30, 45), ("1", 45, 50)]
    >>> range_union(ranges)
    21
    """
    ranges.sort()

    total_len = 0
    cur_chr, cur_left, cur_right = ranges[0] # left-most range
    for r in ranges:
        # open a new range if left > cur_right or chr != cur_chr
        if r[1] > cur_right or r[0] != cur_chr:
            total_len += cur_right - cur_left + 1
            cur_chr, cur_left, cur_right = r
        else:
            # update cur_right
            cur_right = max(r[2], cur_right)

    # the last one
    total_len += cur_right - cur_left + 1

    return total_len


def _make_endpoints(ranges):
    endpoints = []

    for i, (seqid, start, end, score, id) in enumerate(ranges):
        endpoints.append((seqid, start, LEFT, i, score))
        endpoints.append((seqid, end, RIGHT, i, score))

    return sorted(endpoints)


def range_chain(ranges):
    """
    Take a list of weighted intervals, return a non-overlapping set with max weight
    We proceed by looking at the each end point (sorted by their relative positions)

    The input are a list of ranges of the form (start, stop, score), output is a
    subset of the non-overlapping ranges that give the highest score, score

    >>> ranges = [Range("1", 0, 9, 22, 0), Range("1", 3, 18, 24, 1), Range("1", 10, 28, 20, 2)]
    >>> range_chain(ranges)
    ([Range(seqid='1', start=0, end=9, score=22, id=0), Range(seqid='1', start=10, end=28, score=20, id=2)], 42)
    >>> ranges = [Range("2", 0, 1, 3, 0), Range("2", 1, 4, 3, 1), Range("3", 5, 7, 3, 2)]
    >>> range_chain(ranges)
    ([Range(seqid='2', start=0, end=1, score=3, id=0), Range(seqid='3', start=5, end=7, score=3, id=2)], 6)
    """
    ranges.sort()

    endpoints = _make_endpoints(ranges)

    # stores the left end index for quick retrieval
    left_index = {}   
    # dynamic programming, each entry [score, from_index, which_chain]
    scores = []

    for i, (seqid, pos, leftright, j, score) in enumerate(endpoints):

        cur_score = [0, -1, -1] if i==0 else scores[-1][:]

        if leftright == LEFT:
            left_index[j] = i

        else: # this is right end of j-th interval
            # update if chaining j-th interval gives a better score
            left_j = left_index[j]
            chain_score = scores[left_j][0] + score 
            if chain_score > cur_score[0]: 
                cur_score = [chain_score, left_j, j]

        scores.append(cur_score)

    chains = []
    score, last, chain_id = scores[-1] # start backtracking
    while last!=-1:
        if chain_id!=-1: 
            chains.append(chain_id)
        _, last, chain_id = scores[last]

    chains.reverse()

    selected = [ranges[x] for x in chains]

    return selected, score


if __name__ == '__main__':
    
    import doctest
    doctest.testmod()

