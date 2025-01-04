#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Implements the max sum segment algorithm, using Kadane's algorithm, see

<http://en.wikipedia.org/wiki/Maximum_subarray_problem>
"""
Infinity = 1e10000


def max_sum(a):
    """
    For an input array a, output the range that gives the largest sum

    >>> max_sum([4, 4, 9, -5, -6, -1, 5, -6, -8, 9])
    (17, 0, 2)
    >>> max_sum([8, -10, 10, -9, -6, 9, -7, -4, -10, -8])
    (10, 2, 2)
    >>> max_sum([10, 1, -10, -8, 6, 10, -10, 6, -3, 10])
    (19, 4, 9)
    """

    max_sum, max_start_index, max_end_index = -Infinity, 0, 0
    current_max_sum = 0
    current_start_index = 0
    for current_end_index, x in enumerate(a):
        current_max_sum += x
        if current_max_sum > max_sum:
            max_sum, max_start_index, max_end_index = (
                current_max_sum,
                current_start_index,
                current_end_index,
            )
        if current_max_sum < 0:
            current_max_sum = 0
            current_start_index = current_end_index + 1

    return max_sum, max_start_index, max_end_index


if __name__ == "__main__":
    import doctest

    doctest.testmod()

    import numpy as np

    A = np.random.random_integers(-10, 10, 10)
    print("max_sum(%s)" % list(A))
    print(max_sum(A))
