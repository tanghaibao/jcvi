#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Some math formula for various calculations
"""

import sys

from math import log, exp

def jukesCantorD(p, L=100):
    """
    >>> jukesCantorD(.1)
    (0.10732563273050497, 0.0011982248520710061)
    >>> jukesCantorD(.7)
    (2.0310376508266565, 0.47249999999999864)
    """
    assert 0 <= p < .75

    rD = 1 - 4./3 * p
    D = -.75 * log(rD)
    varD = p*(1-p)/(rD**2*L)

    return D, varD


def jukesCantorP(D):
    """
    >>> jukesCantorP(.1)
    0.093620010717789387
    >>> jukesCantorP(2)
    0.6978874115828988
    """
    rD = exp(-4./3 * D)
    p = .75 * (1-rD)
    return p


if __name__ == '__main__':
    import doctest
    doctest.testmod()
