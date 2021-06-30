#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_tsp.py
# algorithms
#
# Created by Haibao Tang on 01/10/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#

import pytest


@pytest.mark.parametrize(
    "edges,directed,expected",
    [
        ([(1, 2), (2, 3), (3, 4), (4, 2), (3, 5)], False, [1, 2, 4, 3, 5]),
        ([(1, 2), (2, 3)], True, [1, 2, 3]),
        (
            [(1, 2, 10), (1, 3, 15), (1, 4, 20), (2, 3, 35), (2, 4, 25), (3, 4, 30)],
            False,
            [3, 1, 2, 4],
        ),
    ],
)
def test_hamiltonian(edges, directed, expected):
    from jcvi.algorithms.tsp import hamiltonian

    solution = hamiltonian(edges, directed, time_limit=1)
    assert solution == expected or solution == list(reversed(expected))
