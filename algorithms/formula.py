#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Some math formula for various calculations
"""

import sys

from math import log, exp
from optparse import OptionParser

from jcvi.utils.cbook import human_size
from jcvi.apps.base import ActionDispatcher, debug
debug()


def jukesCantorD(p, L=100):
    """
    >>> jukesCantorD(.1)
    (0.10732563273050497, 0.0011982248520710061)
    >>> jukesCantorD(.7)
    (2.0310376508266565, 0.47249999999999864)
    """
    assert 0 <= p < .75

    rD = 1 - 4. / 3 * p
    D = -.75 * log(rD)
    varD = p * (1 - p) / (rD ** 2 * L)

    return D, varD


def jukesCantorP(D):
    """
    >>> jukesCantorP(.1)
    0.093620010717789387
    >>> jukesCantorP(2)
    0.6978874115828988
    """
    rD = exp(-4. / 3 * D)
    p = .75 * (1 - rD)
    return p


def main():

    actions = (
        ('velvet', 'calculate velvet memory requirement'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def velvet(args):
    """
    %prog velvet readsize genomesize numreads K

    Calculate velvet memory requirement.
    <http://seqanswers.com/forums/showthread.php?t=2101>

    Ram required for velvetg = -109635 + 18977*ReadSize + 86326*GenomeSize +
    233353*NumReads - 51092*K

    Read size is in bases.
    Genome size is in millions of bases (Mb)
    Number of reads is in millions
    K is the kmer hash value used in velveth
    """
    p = OptionParser(velvet.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    readsize, genomesize, numreads, K = [int(x) for x in args]
    ram = -109635 + 18977 * readsize + 86326 * genomesize + \
            233353 * numreads - 51092 * K
    print >> sys.stderr, "ReadSize: {0}".format(readsize)
    print >> sys.stderr, "GenomeSize: {0}Mb".format(genomesize)
    print >> sys.stderr, "NumReads: {0}M".format(numreads)
    print >> sys.stderr, "K: {0}".format(K)

    ram = human_size(ram * 1000, a_kilobyte_is_1024_bytes=True)
    print >> sys.stderr, "RAM usage: {0} (MAXKMERLENGTH=31)".format(ram)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    main()
