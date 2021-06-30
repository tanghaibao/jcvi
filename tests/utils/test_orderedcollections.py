#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_orderedcollections.py
# utils
#
# Created by Haibao Tang on 03/20/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#

import pytest

from jcvi.utils.orderedcollections import DefaultOrderedDict, parse_qs


@pytest.mark.parametrize(
    "querystr,expected",
    [
        ("a=3;b=4,c=5", DefaultOrderedDict(a=["3"], b=["4,c=5"])),
        ("tt=bb;zz=4,yy=123", DefaultOrderedDict(tt=["bb"], zz=["4,yy=123"])),
    ],
)
def test_parse_qs(querystr, expected):
    attributes = parse_qs(querystr)
    assert attributes == expected
