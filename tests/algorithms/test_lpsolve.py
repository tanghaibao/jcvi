#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_lpsolve.py
# algorithms
#
# Created by Haibao Tang on 01/10/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#

import pytest


def test_has_ortools():
    from jcvi.algorithms.lpsolve import has_ortools

    assert has_ortools()


@pytest.mark.parametrize(
    "constraint_coeffs,bounds,obj_coeffs,num_vars,num_constraints,expected",
    [
        ([{1: 1, 2: 1}], [1], [5, 3, 2], 3, 1, [0, 1]),
        ([{0: 5, 1: 7, 2: 4, 3: 3}], [14], [8, 11, 6, 4], 4, 1, [1, 2, 3]),
    ],
)
def test_solver(
    constraint_coeffs, bounds, obj_coeffs, num_vars, num_constraints, expected
):
    """
    Maximize
     5 x1 + 3 x2 + 2 x3
    Subject to
     x2 + x3 <= 1
    Binary
     x1
     x2
     x3
    End
    """
    from jcvi.algorithms.lpsolve import MIPDataModel

    data = MIPDataModel(
        constraint_coeffs, bounds, obj_coeffs, num_vars, num_constraints
    )
    assert data.solve() == expected
