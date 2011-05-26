#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Matrix related subroutines
"""

import sys

import numpy as np

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()


is_symmetric = lambda M: (M.T == M).all()


def moving_average(a, window=10):
    kernel = np.repeat(1., window) / window
    return np.convolve(a, kernel)[(window-1):]


def determine_signs(nodes, edges):
    """
    Construct the orientation matrix for the pairs on N molecules.

    >>> determine_signs(['A','B','C'],[('A','B','+'),('A','C','-'),('B','C','-')])
    array([ 1,  1, -1])
    """
    N = len(nodes)
    M = np.zeros((N, N), dtype=int)
    for a, b, direction in edges:
        ia, ib = nodes.index(a), nodes.index(b)
        M[ia, ib] = 1 if direction == '+' else -1

    M = symmetrize(M)

    return get_signs(M, validate=False)


def symmetrize(M):
    """
    If M only has a triangle filled with values, all the rest are zeroes,
    this function will copy stuff to the other triangle
    """
    return M + M.T - np.diag(M.diagonal())


def get_signs(M, validate=True):
    """
    Given a numpy array M that contains pairwise orientations, find the largest
    eigenvalue and associated eigenvector and return the signs for the
    eigenvector. This should correspond to the original orientations for the
    individual molecule. In the first example below, let's say 3 molecules A, B
    and C, A-B:same direction, A-C:opposite direction, B-C:opposite
    direction. The final solution is to flip C.

    >>> M = np.array([[0,1,-1],[1,0,-1],[-1,-1,0]])
    >>> get_signs(M)
    array([ 1,  1, -1])
    >>> M = np.array([[0,1,-1],[1,0,0],[-1,0,0]])
    >>> get_signs(M)
    array([ 1,  1, -1])
    """
    # Is this a symmetric matrix?
    assert is_symmetric(M), "the matrix is not symmetric:\n{0}".format(str(M))
    N, x = M.shape

    # eigh() works on symmetric matrix (Hermitian)
    w, v = np.linalg.eigh(M)
    m = np.argmax(w)
    mv = v[:, m]

    sign_array = np.array(np.sign(mv), dtype=int)

    # it does not really matter, but we prefer as few flippings as possible
    if np.sum(sign_array) < 0:
        sign_array = -sign_array

    if validate:
        diag = np.matrix(np.eye(N, dtype=int) * sign_array)
        final = diag * M * diag
        # The final result should have all pairwise in the same direction
        assert (final >= 0).all(), \
                "result check fails:\n{0}".format(final)

    return sign_array


if __name__ == '__main__':
    import doctest
    doctest.testmod()
