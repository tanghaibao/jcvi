#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_gff.py
# formats
#
# Created by Haibao Tang on 03/20/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#

import pytest

from jcvi.formats.gff import DefaultOrderedDict, make_attributes


@pytest.mark.parametrize(
    "attr,gff3,expected",
    [
        (
            "ID=cds00002;Parent=mRNA00002;",
            True,
            DefaultOrderedDict(ID=["cds00002"], Parent=["mRNA00002"]),
        ),
        (
            'Gene 22240.t000374; Note "Carbonic anhydrase"',
            False,
            DefaultOrderedDict(Gene=["22240.t000374"], Note=["Carbonic anhydrase"]),
        ),
    ],
)
def test_make_attributes(attr, gff3, expected):
    features = make_attributes(attr, gff3)
    assert features == expected
