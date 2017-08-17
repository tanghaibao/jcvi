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


ctypedef np.int_t INT
cdef int BB = 20
cdef int *GR = \
     [    842,    1364,    2206,    3571,    5777,    9349,   15126,
        24476,   39602,   64079,  103681,  167761,  271442,  439204,
       710646, 1149851, 1860497, 3010349, 4870846, 7881196]


@cython.boundscheck(False)  # Turn off bounds-checking
@cython.wraparound(False)   # Turn off negative index wrapping
def score_evaluate(array.array tour,
                   np.ndarray[INT, ndim=1] tour_sizes=None,
                   np.ndarray[INT, ndim=2] tour_M=None):

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


@cython.boundscheck(False)  # Turn off bounds-checking
@cython.wraparound(False)   # Turn off negative index wrapping
def score_evaluate_oriented(array.array tour,
                            np.ndarray[INT, ndim=1] tour_sizes=None,
                            np.ndarray[INT, ndim=3] tour_P=None):
    cdef np.ndarray sizes_oo = \
                np.array([tour_sizes[x] for x in tour])
    cdef np.ndarray sizes_cum = np.cumsum(sizes_oo)

    cdef double s = 0.0
    cdef int size = len(tour)
    cdef int a, b, c, ia, ib, ic
    for ia in range(size):
        a = tour[ia]
        for ib in range(ia + 1, size):
            b = tour[ib]
            dist = sizes_cum[ib - 1] - sizes_cum[ia]
            if dist > 1e7:
                break
            for ic in range(BB):
                c = tour_P[a, b, ic]
                s += c / (GR[ic] + dist)
    return s,
