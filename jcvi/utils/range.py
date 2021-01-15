#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This script implements algorithm for finding intersecting rectangles,
both on the 2D dotplot and 1D-projection

`range_chain` implements the exon-chain algorithm
"""
from __future__ import print_function

import sys

from collections import namedtuple, defaultdict
from itertools import groupby
from more_itertools import pairwise


LEFT, RIGHT = 0, 1
Range = namedtuple("Range", "seqid start end score id")


def range_parse(s):
    """
    >>> range_parse("chr1:1000-1")
    Range(seqid='chr1', start=1, end=1000, score=0, id=0)
    """
    chr, se = s.split(":")
    start, end = se.split("-")
    start, end = int(start), int(end)
    if start > end:
        start, end = end, start

    return Range(chr, start, end, 0, 0)


def range_intersect(a, b, extend=0):
    """
    Returns the intersection between two reanges.

    >>> range_intersect((30, 45), (55, 65))
    >>> range_intersect((48, 65), (45, 55))
    [48, 55]
    """
    a_min, a_max = a
    if a_min > a_max:
        a_min, a_max = a_max, a_min
    b_min, b_max = b
    if b_min > b_max:
        b_min, b_max = b_max, b_min

    if a_max + extend < b_min or b_max + extend < a_min:
        return None
    i_min = max(a_min, b_min)
    i_max = min(a_max, b_max)
    if i_min > i_max + extend:
        return None

    return [i_min, i_max]


def ranges_intersect(rset):
    """
    Recursively calls the range_intersect() - pairwise version.

    >>> ranges_intersect([(48, 65), (45, 55), (50, 56)])
    [50, 55]
    """
    if not rset:
        return None

    a = rset[0]
    for b in rset[1:]:
        if not a:
            return None
        a = range_intersect(a, b)

    return a


def range_overlap(a, b, ratio=False):
    """
    Returns whether two ranges overlap. Set percentage=True returns overlap
    ratio over the shorter range of the two.

    >>> range_overlap(("1", 30, 45), ("1", 41, 55))
    5
    >>> range_overlap(("1", 21, 45), ("1", 41, 75), ratio=True)
    0.2
    >>> range_overlap(("1", 30, 45), ("1", 15, 55))
    16
    >>> range_overlap(("1", 30, 45), ("1", 15, 55), ratio=True)
    1.0
    >>> range_overlap(("1", 30, 45), ("1", 57, 68))
    0
    >>> range_overlap(("1", 30, 45), ("2", 42, 55))
    0
    >>> range_overlap(("1", 30, 45), ("2", 42, 55), ratio=True)
    0.0
    """
    a_chr, a_min, a_max = a
    b_chr, b_min, b_max = b
    a_min, a_max = sorted((a_min, a_max))
    b_min, b_max = sorted((b_min, b_max))
    shorter = min((a_max - a_min), (b_max - b_min)) + 1
    # must be on the same chromosome
    if a_chr != b_chr:
        ov = 0
    else:
        ov = min(shorter, (a_max - b_min + 1), (b_max - a_min + 1))
        ov = max(ov, 0)
    if ratio:
        ov /= float(shorter)
    return ov


def range_distance(a, b, distmode="ss"):
    """
    Returns the distance between two ranges.

    distmode is ss, se, es, ee and sets the place on read one and two to
          measure the distance (s = start, e = end)

    >>> range_distance(("1", 30, 45, '+'), ("1", 45, 55, '+'))
    (26, '++')
    >>> range_distance(("1", 30, 45, '-'), ("1", 57, 68, '-'))
    (39, '--')
    >>> range_distance(("1", 30, 42, '-'), ("1", 45, 55, '+'))
    (26, '-+')
    >>> range_distance(("1", 30, 42, '+'), ("1", 45, 55, '-'), distmode='ee')
    (2, '+-')
    """
    assert distmode in ("ss", "ee")

    a_chr, a_min, a_max, a_strand = a
    b_chr, b_min, b_max, b_strand = b
    # must be on the same chromosome
    if a_chr != b_chr:
        dist = -1
    # elif range_overlap(a[:3], b[:3]):
    #    dist = 0
    else:
        # If the two ranges do not overlap, check stranded-ness and distance
        if a_min > b_min:
            a_min, b_min = b_min, a_min
            a_max, b_max = b_max, a_max
            a_strand, b_strand = b_strand, a_strand

        if distmode == "ss":
            dist = b_max - a_min + 1
        elif distmode == "ee":
            dist = b_min - a_max - 1

    orientation = a_strand + b_strand

    return dist, orientation


def range_minmax(ranges):
    """
    Returns the span of a collection of ranges where start is the smallest of
    all starts, and end is the largest of all ends.

    >>> ranges = [(30, 45), (40, 50), (10, 100)]
    >>> range_minmax(ranges)
    (10, 100)
    """
    rmin = min(ranges)[0]
    rmax = max(ranges, key=lambda x: x[1])[1]
    return rmin, rmax


def range_closest(ranges, b, left=True):
    """
    Returns the range that's closest to the given position. Notice that the
    behavior is to return ONE closest range to the left end (if left is True).
    This is a SLOW method.

    >>> ranges = [("1", 30, 40), ("1", 33, 35), ("1", 10, 20)]
    >>> b = ("1", 22, 25)
    >>> range_closest(ranges, b)
    ('1', 10, 20)
    >>> range_closest(ranges, b, left=False)
    ('1', 33, 35)
    >>> b = ("1", 2, 5)
    >>> range_closest(ranges, b)
    """
    from jcvi.utils.orderedcollections import SortedCollection

    key = (lambda x: x) if left else (lambda x: (x[0], x[2], x[1]))
    rr = SortedCollection(ranges, key=key)
    try:
        if left:
            s = rr.find_le(b)
            assert key(s) <= key(b), (s, b)
        else:
            s = rr.find_ge(b)
            assert key(s) >= key(b), (s, b)
    except ValueError:
        s = None

    return s


def range_interleave(ranges, sizes={}, empty=False):
    """
    Returns the ranges in between the given ranges.

    >>> ranges = [("1", 30, 40), ("1", 45, 50), ("1", 10, 30)]
    >>> range_interleave(ranges)
    [('1', 41, 44)]
    >>> ranges = [("1", 30, 40), ("1", 42, 50)]
    >>> range_interleave(ranges)
    [('1', 41, 41)]
    >>> range_interleave(ranges, sizes={"1": 70})
    [('1', 1, 29), ('1', 41, 41), ('1', 51, 70)]
    """
    ranges = range_merge(ranges)
    interleaved_ranges = []

    for ch, cranges in groupby(ranges, key=lambda x: x[0]):
        cranges = list(cranges)
        size = sizes.get(ch, None)
        if size:
            ch, astart, aend = cranges[0]
            if astart > 1:
                interleaved_ranges.append((ch, 1, astart - 1))
            elif empty:
                interleaved_ranges.append(None)

        for a, b in pairwise(cranges):
            ch, astart, aend = a
            ch, bstart, bend = b
            istart, iend = aend + 1, bstart - 1
            if istart <= iend:
                interleaved_ranges.append((ch, istart, iend))
            elif empty:
                interleaved_ranges.append(None)

        if size:
            ch, astart, aend = cranges[-1]
            if aend < size:
                interleaved_ranges.append((ch, aend + 1, size))
            elif empty:
                interleaved_ranges.append(None)

    return interleaved_ranges


def range_merge(ranges, dist=0):
    """
    Returns merged range. Similar to range_union, except this returns
    new ranges.

    >>> ranges = [("1", 30, 45), ("1", 40, 50), ("1", 10, 50)]
    >>> range_merge(ranges)
    [('1', 10, 50)]
    >>> ranges = [("1", 30, 40), ("1", 45, 50)]
    >>> range_merge(ranges)
    [('1', 30, 40), ('1', 45, 50)]
    >>> ranges = [("1", 30, 40), ("1", 45, 50)]
    >>> range_merge(ranges, dist=5)
    [('1', 30, 50)]
    """
    if not ranges:
        return []

    ranges.sort()

    cur_range = list(ranges[0])
    merged_ranges = []
    for r in ranges[1:]:
        # open new range if start > cur_end or seqid != cur_seqid
        if r[1] - cur_range[2] > dist or r[0] != cur_range[0]:
            merged_ranges.append(tuple(cur_range))
            cur_range = list(r)
        else:
            cur_range[2] = max(cur_range[2], r[2])
    merged_ranges.append(tuple(cur_range))

    return merged_ranges


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
    >>> range_union([])
    0
    """
    if not ranges:
        return 0

    ranges.sort()

    total_len = 0
    cur_chr, cur_left, cur_right = ranges[0]  # left-most range
    for r in ranges:
        # open new range if left > cur_right or chr != cur_chr
        if r[1] > cur_right or r[0] != cur_chr:
            total_len += cur_right - cur_left + 1
            cur_chr, cur_left, cur_right = r
        else:
            # update cur_right
            cur_right = max(r[2], cur_right)

    # the last one
    total_len += cur_right - cur_left + 1

    return total_len


def range_span(ranges):
    """
    Returns the total span between the left most range to the right most range.

    >>> ranges = [("1", 30, 45), ("1", 40, 50), ("1", 10, 50)]
    >>> range_span(ranges)
    41
    >>> ranges = [("1", 30, 45), ("2", 40, 50)]
    >>> range_span(ranges)
    27
    >>> ranges = [("1", 30, 45), ("1", 45, 50)]
    >>> range_span(ranges)
    21
    >>> range_span([])
    0
    """
    if not ranges:
        return 0

    ranges.sort()
    ans = 0
    for seq, lt in groupby(ranges, key=lambda x: x[0]):
        lt = list(lt)
        ans += max(max(lt)[1:]) - min(min(lt)[1:]) + 1
    return ans


def _make_endpoints(ranges):
    assert ranges, "Ranges cannot be empty"
    endpoints = []

    for i, (seqid, start, end, score, id) in enumerate(ranges):
        endpoints.append((seqid, start, LEFT, i, score))
        endpoints.append((seqid, end, RIGHT, i, score))

    return sorted(endpoints)


def range_piles(ranges):
    """
    Return piles of intervals that overlap. The piles are only interrupted by
    regions of zero coverage.

    >>> ranges = [Range("2", 0, 1, 3, 0), Range("2", 1, 4, 3, 1), Range("3", 5, 7, 3, 2)]
    >>> list(range_piles(ranges))
    [[0, 1], [2]]
    """
    endpoints = _make_endpoints(ranges)

    for seqid, ends in groupby(endpoints, lambda x: x[0]):
        active = []
        depth = 0
        for seqid, pos, leftright, i, score in ends:
            if leftright == LEFT:
                active.append(i)
                depth += 1
            else:
                depth -= 1

            if depth == 0 and active:
                yield active
                active = []


def range_conflict(ranges, depth=1):
    """
    Find intervals that are overlapping in 1-dimension.
    Return groups of block IDs that are in conflict.

    >>> ranges = [Range("2", 0, 1, 3, 0), Range("2", 1, 4, 3, 1), Range("3", 5, 7, 3, 2)]
    >>> list(range_conflict(ranges))
    [(0, 1)]
    """
    overlap = set()
    active = set()
    endpoints = _make_endpoints(ranges)

    for seqid, ends in groupby(endpoints, lambda x: x[0]):
        active.clear()
        for seqid, pos, leftright, i, score in ends:
            if leftright == LEFT:
                active.add(i)
            else:
                active.remove(i)

            if len(active) > depth:
                overlap.add(tuple(sorted(active)))

    for ov in overlap:
        yield ov


def range_chain(ranges):
    """
    Take list of weighted intervals, find non-overlapping set with max weight.
    We proceed with each end point (sorted by their relative positions).

    The input are a list of ranges of the form (start, stop, score), output is
    subset of the non-overlapping ranges that give the highest score, score

    >>> ranges = [Range("1", 0, 9, 22, 0), Range("1", 3, 18, 24, 1), Range("1", 10, 28, 20, 2)]
    >>> range_chain(ranges)
    ([Range(seqid='1', start=0, end=9, score=22, id=0), Range(seqid='1', start=10, end=28, score=20, id=2)], 42)
    >>> ranges = [Range("2", 0, 1, 3, 0), Range("2", 1, 4, 3, 1), Range("3", 5, 7, 3, 2)]
    >>> range_chain(ranges)
    ([Range(seqid='2', start=0, end=1, score=3, id=0), Range(seqid='3', start=5, end=7, score=3, id=2)], 6)
    """
    endpoints = _make_endpoints(ranges)

    # stores the left end index for quick retrieval
    left_index = {}
    # dynamic programming, each entry [score, from_index, which_chain]
    scores = []

    for i, (seqid, pos, leftright, j, score) in enumerate(endpoints):

        cur_score = [0, -1, -1] if i == 0 else scores[-1][:]

        if leftright == LEFT:
            left_index[j] = i

        else:  # this is right end of j-th interval
            # update if chaining j-th interval gives a better score
            left_j = left_index[j]
            chain_score = scores[left_j][0] + score
            if chain_score > cur_score[0]:
                cur_score = [chain_score, left_j, j]

        scores.append(cur_score)

    chains = []
    score, last, chain_id = scores[-1]  # start backtracking
    while last != -1:
        if chain_id != -1:
            chains.append(chain_id)
        _, last, chain_id = scores[last]

    chains.reverse()

    selected = [ranges[x] for x in chains]

    return selected, score


def ranges_depth(ranges, sizes, verbose=True):
    """
    Allow triple (seqid, start, end) rather than just tuple (start, end)
    """
    ranges.sort()
    for seqid, rrs in groupby(ranges, key=lambda x: x[0]):
        rrs = [(a, b) for (s, a, b) in rrs]
        size = sizes[seqid]
        ds, depthdetails = range_depth(rrs, size, verbose=verbose)
        depthdetails = [(seqid, s, e, d) for s, e, d in depthdetails]
        yield depthdetails


def range_depth(ranges, size, verbose=True):
    """
    Overlay ranges on [start, end], and summarize the ploidy of the intervals.
    """
    from jcvi.utils.cbook import percentage

    # Make endpoints
    endpoints = []
    for a, b in ranges:
        endpoints.append((a, LEFT))
        endpoints.append((b, RIGHT))
    endpoints.sort()
    vstart, vend = min(endpoints)[0], max(endpoints)[0]

    assert 0 <= vstart < size
    assert 0 <= vend < size

    depth = 0
    depthstore = defaultdict(int)
    depthstore[depth] += vstart
    depthdetails = [(0, vstart, depth)]

    for (a, atag), (b, btag) in pairwise(endpoints):
        if atag == LEFT:
            depth += 1
        elif atag == RIGHT:
            depth -= 1
        depthstore[depth] += b - a
        depthdetails.append((a, b, depth))

    assert btag == RIGHT
    depth -= 1

    assert depth == 0
    depthstore[depth] += size - vend
    depthdetails.append((vend, size, depth))

    assert sum(depthstore.values()) == size
    if verbose:
        for depth, count in sorted(depthstore.items()):
            print(
                "Depth {0}: {1}".format(depth, percentage(count, size)), file=sys.stderr
            )

    return depthstore, depthdetails


if __name__ == "__main__":

    import doctest

    doctest.testmod()
