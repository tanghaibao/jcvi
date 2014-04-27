#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Multiple TSPs problem are cases where several salesmen travel together. As they
travel between cities, they each incur different cost. In addition, they visit
different sets cities, but are constrained by a common visiting order. For
illustration:

            A B C D E F G H
Salesman 1  x   x   x x   x
Salesman 2    x x x     x x

The goal is to minimize total combined cost for all salesmen. This is sometimes
also called "synchronized TSP". In the above example, final cost is:

Salesman 1 cost: AC + CE + EF + FH
Salesman 2 cost: BC + CD + DG + GH
"""

from jcvi.utils.iter import pairwise
from jcvi.algorithms.lpsolve import lpsolve_tsp_tour
from jcvi.algorithms.tsp import concorde_tour, make_data


def evaluate(tour, M):
    score = 0
    for ia, ib in pairwise(tour):
        score += M[ia, ib]
    return score


def main():
    POINTS = 20
    x, y, M = make_data(POINTS)
    ltour = lpsolve_tsp_tour(POINTS, M)
    print ltour, evaluate(ltour, M)

    ctour = concorde_tour(POINTS, M)
    print ctour, evaluate(ctour, M)


if __name__ == '__main__':
    main()
