#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_dotplot.py
# graphics
#
# Created by Haibao Tang on 06/19/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#


def test_downsample():
    from random import random
    from jcvi.graphics.dotplot import downsample

    lots_of_data = [random() for _ in range(20000)]
    assert len(lots_of_data) == 20000
    assert len(downsample(lots_of_data)) == 10000
    assert len(downsample(lots_of_data, 5000)) == 5000

    few_data = [random() for _ in range(10)]
    assert len(few_data) == 10
    assert len(downsample(few_data)) == 10
    assert len(downsample(few_data, 1)) == 1
