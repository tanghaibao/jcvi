#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Cythonized version of score_evaluate() in hic.py.
"""

from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
from cpython cimport array
import array


ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)  # Turn off bounds-checking
@cython.wraparound(False)   # Turn off negative index wrapping for entire
def score_evaluate(array.array tour,
                   np.ndarray[DTYPE_t, ndim=1] tour_sizes=None,
                   np.ndarray[DTYPE_t, ndim=2] tour_M=None):

    cdef int x
    cdef np.ndarray sizes_oo = \
                np.array([tour_sizes[x] for x in tour])
    cdef np.ndarray sizes_cum = np.cumsum(sizes_oo) - sizes_oo // 2

    cdef double s = 0.0
    cdef int size = len(tour)
    cdef int a, b, links
    cdef int ia, ib
    cdef int dist
    for ia in range(size):
        a = tour[ia]
        for ib in range(ia + 1, size):
            b = tour[ib]
            links = tour_M[a, b]
            dist = sizes_cum[ib] - sizes_cum[ia]
            if dist > 1e7:
                break
            s += links / dist
    return s,
