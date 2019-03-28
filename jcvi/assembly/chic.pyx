#cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

"""
Cythonized version of score_evaluate() in hic.py.

Support three versions with different objective functions:
- score_evaluate_M: distance is defined as the distance between mid-points
  between contigs. Maximize Sum(n_links / distance).
- score_evaluate_P: distance is defined as the sizes of interleaving contigs
  plus the harmonic mean of all link distances. Maximize Sum(n_links / distance).
- score_evaluate_Q: distance is defined as the sizes of interleaving contigs
  plus the actual link distances. Maximize Sum(1 / distance) for all links.
  For performance consideration, we actually use a histogram to approximate
  all link distances. See golden_array() in hic for details.
"""

from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
from cpython cimport array
import array


ctypedef np.int_t INT
DEF LIMIT = 10000000
DEF BB = 12
cdef int *GR = \
     [   5778,    9349,   15127,   24476,
        39603,   64079,  103682,  167761,
       271443,  439204,  710647, 1149851]


def score_evaluate_M(array.array[int] tour,
                     np.ndarray[INT, ndim=1] tour_sizes=None,
                     np.ndarray[INT, ndim=2] tour_M=None):
    cdef np.ndarray[INT, ndim=1] sizes_oo = tour_sizes[tour]
    cdef np.ndarray[INT, ndim=1] sizes_cum = np.cumsum(sizes_oo) - sizes_oo // 2

    cdef double s = 0.0
    cdef int size = len(tour)
    cdef int a, b, ia, ib
    cdef int links
    cdef double dist
    for ia in range(size):
        a = tour[ia]
        for ib in range(ia + 1, size):
            b = tour[ib]
            links = tour_M[a, b]
            if links == 0:
                continue
            dist = sizes_cum[ib] - sizes_cum[ia]
            if dist > LIMIT:
                break
            s += links / dist
    return s,


def score_evaluate_P(array.array[int] tour,
                     np.ndarray[INT, ndim=1] tour_sizes=None,
                     np.ndarray[INT, ndim=3] tour_P=None):
    cdef np.ndarray[INT, ndim=1] sizes_oo = tour_sizes[tour]
    cdef np.ndarray[INT, ndim=1] sizes_cum = np.cumsum(sizes_oo)

    cdef double s = 0.0
    cdef int size = len(tour)
    cdef int a, b, c, ia, ib
    cdef double dist
    for ia in range(size):
        a = tour[ia]
        for ib in range(ia + 1, size):
            b = tour[ib]
            dist = sizes_cum[ib - 1] - sizes_cum[ia]
            if dist > LIMIT:
                break
            c = tour_P[a, b, 0]
            if c == 0:
                continue
            s += c / (tour_P[a, b, 1] + dist)
    return s,


def score_evaluate_Q(array.array[int] tour,
                     np.ndarray[INT, ndim=1] tour_sizes=None,
                     np.ndarray[INT, ndim=3] tour_Q=None):
    cdef np.ndarray[INT, ndim=1] sizes_oo = tour_sizes[tour]
    cdef np.ndarray[INT, ndim=1] sizes_cum = np.cumsum(sizes_oo)

    cdef double s = 0.0
    cdef int size = len(tour)
    cdef int a, b, c, ia, ib, ic
    cdef double dist
    for ia in range(size):
        a = tour[ia]
        for ib in range(ia + 1, size):
            b = tour[ib]
            if tour_Q[a, b, 0] == -1:
                continue
            dist = sizes_cum[ib - 1] - sizes_cum[ia]
            if dist > LIMIT:
                break
            for ic in range(BB):
                c = tour_Q[a, b, ic]
                s += c / (GR[ic] + dist)
    return s,
