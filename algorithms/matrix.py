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


is_symmetric = lambda M: (M.transpose()==M).all()

def get_signs(M, validate=True):
    """
    Given a numpy array M that contains pairwise orientations, find the largest
    eigenvalue and associated eigenvector and return the signs for the
    eigenvector. This should correspond to the original orientations for the
    individual molecule. In the first example below, let's say 3 molecules A, B
    and C, A-B:opposite direction, A-C:opposite direction, B-C: positive
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

    w, v = np.linalg.eig(M)
    m = np.argmax(w)
    mv = v[:, m]
    sign_array = np.array(np.sign(mv), dtype=int)

    if validate:
        diag = np.matrix(np.eye(3, dtype=int) * sign_array)
        final = diag * M * diag
        # The final result should have all pairwise in the same direction
        assert (final>=0).all(), \
                "result check fails:\n{0}\n{1}".format(final, target)

    return sign_array


if __name__ == '__main__':
    import doctest
    doctest.testmod()
