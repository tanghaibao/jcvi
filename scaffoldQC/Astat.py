#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Discriminator A-statistics:

If n reads are uniform sample of the genome of length G, 
we expect k = n * delta / G to start in a region of length delta

Use poisson distribution:
A(delta, k) = ln(prob(1-copy) / prob(2-copies)) = n * delta / G - k * ln2
"""

from math import log

ln2 = log(2)


def Astat(delta, k, G, n):
    """
    delta: contig size
    k: reads mapped in contig
    G: total genome size
    n: total reads mapped to genome
    """
    return n * delta * 1. / G - k * ln2

