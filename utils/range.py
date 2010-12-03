#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This script implements algorithm for finding intersecting rectangles, 
both on the 2D dotplot and 1D-projection
"""

from grouper import Grouper

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


def range_union(ranges):
    """
    Returns total size of ranges, expect range as (chr, left, right)
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


def range_chain(ranges):
    """
    ranges in the form of (start, stop, score), get the best scoring chain

    >>> ranges = [(0, 9, 22), (3, 18, 24), (10, 28, 20), (20, 39, 26)]
    >>> range_chain(ranges)
    50
    """
    ranges.sort()

    endpoints = []
    LEFT, RIGHT = 0, 1

    for i, (start, stop, score) in enumerate(ranges):
        endpoints.append((start, i, score, LEFT))
        endpoints.append((stop, i, score, RIGHT))

    endpoints.sort()

    left_index = {} 
    scores = [0]
    cur_score = 0

    for i, (pos, j, score, leftright) in enumerate(endpoints):

        if leftright == LEFT:
            left_index[j] = i + 1

        else:
            left_j = left_index[j]
            chain_score = scores[left_j] + score
            if chain_score > cur_score:
                cur_score = chain_score

        scores.append(cur_score)

    return scores[-1]


if __name__ == '__main__':
    
    import doctest
    doctest.testmod()

