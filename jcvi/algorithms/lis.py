#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Longest increasing subsequence, code stolen from internet (thanks)
# http://wordaligned.org/articles/patience-sort
from __future__ import print_function
import bisect

# We want a maximum function which accepts a default value
from functools import partial, reduce

maximum = partial(reduce, max)


def patience_sort(xs):
    """Patience sort an iterable, xs.

    This function generates a series of pairs (x, pile), where "pile"
    is the 0-based index of the pile "x" should be placed on top of.
    Elements of "xs" must be less-than comparable.
    """
    pile_tops = list()
    for x in xs:
        pile = bisect.bisect_left(pile_tops, x)
        if pile == len(pile_tops):
            pile_tops.append(x)
        else:
            pile_tops[pile] = x
        yield x, pile


def longest_monotonic_subseq_length(xs):
    """Return the length of the longest monotonic subsequence of xs, second
    return value is the difference between increasing and decreasing lengths.

    >>> longest_monotonic_subseq_length((4, 5, 1, 2, 3))
    (3, 1)
    >>> longest_monotonic_subseq_length((1, 2, 3, 5, 4))
    (4, 2)
    >>> longest_monotonic_subseq_length((1, 2, 1))
    (2, 0)
    """
    li = longest_increasing_subseq_length(xs)
    ld = longest_decreasing_subseq_length(xs)
    return max(li, ld), li - ld


def longest_increasing_subseq_length(xs):
    """Return the length of the longest increasing subsequence of xs.

    >>> longest_increasing_subseq_length(range(3))
    3
    >>> longest_increasing_subseq_length([3, 1, 2, 0])
    2
    """
    return 1 + maximum((pile for x, pile in patience_sort(xs)), -1)


def longest_decreasing_subseq_length(xs):
    return longest_increasing_subseq_length(reversed(xs))


def longest_monotonic_subseq_length_loose(xs):
    li = longest_increasing_subseq_length_loose(xs)
    ld = longest_decreasing_subseq_length_loose(xs)
    return max(li, ld), li - ld


def longest_increasing_subseq_length_loose(xs):
    xs = [(x, i) for (i, x) in enumerate(xs)]
    return longest_increasing_subseq_length(xs)


def longest_decreasing_subseq_length_loose(xs):
    xs = [(x, -i) for (i, x) in enumerate(xs)]
    return longest_decreasing_subseq_length(xs)


def longest_increasing_subsequence(xs):
    """Return a longest increasing subsequence of xs.

    (Note that there may be more than one such subsequence.)
    >>> longest_increasing_subsequence(range(3))
    [0, 1, 2]
    >>> longest_increasing_subsequence([3, 1, 2, 0])
    [1, 2]
    """
    # Patience sort xs, stacking (x, prev_ix) pairs on the piles.
    # Prev_ix indexes the element at the top of the previous pile,
    # which has a lower x value than the current x value.
    piles = [[]]  # Create a dummy pile 0
    for x, p in patience_sort(xs):
        if p + 1 == len(piles):
            piles.append([])
        # backlink to the top of the previous pile
        piles[p + 1].append((x, len(piles[p]) - 1))
    # Backtrack to find a longest increasing subsequence
    npiles = len(piles) - 1
    prev = 0
    lis = list()
    for pile in range(npiles, 0, -1):
        x, prev = piles[pile][prev]
        lis.append(x)
    lis.reverse()
    return lis


def longest_decreasing_subsequence(xs):
    """
    Wrapper that calls longest_increasing_subsequence
    >>> longest_decreasing_subsequence([23, 19, 97, 16, 37, 44, 88, 77, 26])
    [97, 88, 77, 26]
    """
    return list(reversed(longest_increasing_subsequence(reversed(xs))))


def longest_monotonic_subsequence(xs):
    lis = longest_increasing_subsequence(xs)
    lds = longest_decreasing_subsequence(xs)
    if len(lis) >= len(lds):
        return lis
    return lds


def longest_monotonic_subsequence_loose(xs):
    lis = longest_increasing_subsequence_loose(xs)
    lds = longest_decreasing_subsequence_loose(xs)
    if len(lis) >= len(lds):
        return lis
    return lds


def longest_increasing_subsequence_loose(xs):
    xs = [(x, i) for (i, x) in enumerate(xs)]
    ll = longest_increasing_subsequence(xs)
    return [x for (x, i) in ll]


def longest_decreasing_subsequence_loose(xs):
    xs = [(x, -i) for (i, x) in enumerate(xs)]
    ll = longest_decreasing_subsequence(xs)
    return [x for (x, i) in ll]


def backtracking(a, L, bestsofar):
    """
    Start with the heaviest weight and emit index
    """
    w, j = max(L.items())
    while j != -1:
        yield j
        w, j = bestsofar[j]


def heaviest_increasing_subsequence(a, debug=False):
    """
    Returns the heaviest increasing subsequence for array a. Elements are (key,
    weight) pairs.

    >>> heaviest_increasing_subsequence([(3, 3), (2, 2), (1, 1), (0, 5)])
    ([(0, 5)], 5)
    """
    # Stores the smallest idx of last element of a subsequence of weight w
    L = {0: -1}
    bestsofar = [(0, -1)] * len(a)  # (best weight, from_idx)
    for i, (key, weight) in enumerate(a):

        for w, j in list(L.items()):
            if j != -1 and a[j][0] >= key:
                continue

            new_weight = w + weight
            if new_weight in L and a[L[new_weight]][0] <= key:
                continue

            L[new_weight] = i
            newbest = (new_weight, j)
            if newbest > bestsofar[i]:
                bestsofar[i] = newbest

        if debug:
            # print (key, weight), L
            print((key, weight), bestsofar)

    tb = reversed(list(backtracking(a, L, bestsofar)))
    return [a[x] for x in tb], max(L.items())[0]


if __name__ == "__main__":
    import doctest

    doctest.testmod()

    import numpy as np

    LENGTH = 20
    A = [np.random.randint(0, 20) for x in range(LENGTH)]
    A = list(A)
    B = list(zip(A, [1] * LENGTH))
    print(A)
    lis = longest_increasing_subsequence(A)
    print("longest increasing:", lis)
    lds = longest_decreasing_subsequence(A)
    print("longest decreasing:", lds)
    lisl = longest_increasing_subsequence_loose(A)
    print("longest increasing loose:", lisl)
    ldsl = longest_decreasing_subsequence_loose(A)
    print("longest decreasing loose:", ldsl)
    # this should be the same as longest_increasing_subsequence
    his, his_dd = heaviest_increasing_subsequence(B)
    hlis, wts = zip(*his)
    print("heaviest increasing (weight 1, compare with lis):", hlis)
    assert len(lis) == len(his)
