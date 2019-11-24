#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_matrix.py
# algorithms
#
# Created by Haibao Tang on 11/24/19
# Copyright © 2019 Haibao Tang. All rights reserved.
#

import pytest
import numpy as np


@pytest.mark.parametrize(
    "input_array,factor,expected",
    [
        (np.arange(16, dtype="int").reshape(4, 4), 2, np.array([[10, 18], [42, 50]])),
        (np.arange(16, dtype="int").reshape(4, 4), 4, np.array([[120]])),
    ],
)
def test_compact(input_array, factor, expected):
    from jcvi.algorithms.matrix import compact

    assert np.allclose(compact(input_array, factor=factor), expected)


@pytest.mark.parametrize(
    "nodes,edges,expected",
    [([0, 1, 2], [(0, 1, 1), (1, 2, 2), (0, 2, 3)], np.array([0, 1, 3]))],
)
def test_determine_positions(nodes, edges, expected):
    from jcvi.algorithms.matrix import determine_positions

    assert np.allclose(determine_positions(nodes, edges), expected)


@pytest.mark.parametrize(
    "nodes,edges,expected",
    [([0, 1, 2], [(0, 1, 1), (0, 2, -1), (1, 2, -1)], np.array([1, 1, -1]))],
)
def test_determine_signs(nodes, edges, expected):
    from jcvi.algorithms.matrix import determine_signs

    assert np.allclose(determine_signs(nodes, edges), expected)


@pytest.mark.parametrize(
    "input_array,expected",
    [
        (np.array([[0, 1, -1], [1, 0, -1], [-1, -1, 0]]), np.array([1, 1, -1])),
        (np.array([[0, 1, -1], [1, 0, 0], [-1, 0, 0]]), np.array([1, 1, -1])),
    ],
)
def test_get_signs(input_array, expected):
    from jcvi.algorithms.matrix import get_signs

    assert np.allclose(get_signs(input_array), expected)


@pytest.mark.parametrize(
    "A,K,L,expected",
    [
        (
            np.array([[1, -1, 0], [0, 1, -1], [1, 0, -1]]),
            np.eye(3, dtype=int),
            np.array([1, 2, 3]),
            np.array([1.0, 3.0]),
        )
    ],
)
def test_spring_system(A, K, L, expected):
    from jcvi.algorithms.matrix import spring_system

    assert np.allclose(spring_system(A, K, L), expected)
