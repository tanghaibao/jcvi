#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_maxsum.py
# algorithms
#
# Created by Haibao Tang on 06/19/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#

import pytest


@pytest.mark.parametrize(
    "input,output",
    [
        ([4, 4, 9, -5, -6, -1, 5, -6, -8, 9], (17, 0, 2)),
        ([8, -10, 10, -9, -6, 9, -7, -4, -10, -8], (10, 2, 2)),
        ([10, 1, -10, -8, 6, 10, -10, 6, -3, 10], (19, 4, 9)),
    ],
)
def test_max_sum(input, output):
    from jcvi.algorithms.maxsum import max_sum

    assert max_sum(input) == output
