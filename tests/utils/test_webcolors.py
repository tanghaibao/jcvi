#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_webcolors.py
# utils
#
# Created by Haibao Tang on 01/09/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#

import pytest


@pytest.mark.parametrize(
    "rgb1,rgb2,output",
    [
        ((0, 0, 0), (0, 0, 0), 0),
        ((0, 0, 0), (0, 0, 1), 0.6322079321995702),
        ((10, 20, 30), (110, 120, 130), 43.07089640618708),
    ],
)
def test_color_diff(rgb1, rgb2, output):
    from jcvi.utils.webcolors import color_diff

    assert color_diff(rgb1, rgb2) == pytest.approx(output)


@pytest.mark.parametrize(
    "input,output",
    [((171, 222, 230), "powderblue"), ((254, 225, 232), "lavenderblush")],
)
def test_closest_color(input, output):
    from jcvi.utils.webcolors import closest_color

    assert closest_color(input) == output
