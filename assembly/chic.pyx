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
cdef int BB = 12
cdef int *GR = \
     [   5778,    9349,   15127,   24476,
        39603,   64079,  103682,  167761,
       271443,  439204,  710647, 1149851]


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
    cdef int a, b, ia, ib
    cdef int links
    cdef int dist
    for ia in range(size):
        a = tour[ia]
        for ib in range(ia + 1, size):
            b = tour[ib]
            links = tour_M[a, b]
            if links == 0:
                continue
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
    cdef int dist
    for ia in range(size):
        a = tour[ia]
        for ib in range(ia + 1, size):
            b = tour[ib]
            if tour_P[a, b, 0] == -1:
                continue
            dist = sizes_cum[ib - 1] - sizes_cum[ia]
            if dist > 1e7:
                break
            for ic in range(BB):
                c = tour_P[a, b, ic]
                s += c / (GR[ic] + dist)
    return s,
