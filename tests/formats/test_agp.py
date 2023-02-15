#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_agp.py
# formats
#
# Created by Haibao Tang on 12/25/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#

import pytest


@pytest.mark.parametrize(
    "line, should_assert",
    [
        ("chr23\t9512061\t10491479\t25\tW\tscaffold_164\t1\t979419\t+", False),
        ("chr23\t10491480\t10491579\t26\tU\t100\tbad_gaptype\tyes", True),
        (
            "chr23\t10491480\t10491579\t26\tU\t100\tcontig\tyes\tproximity_ligation",
            False,
        ),
        ("chr23\t10491480\t10491579\t26\tU\t100\tcontig\tyes\tbad_linkage", True),
    ],
)
def test_agp_validate(line: str, should_assert: bool):
    from jcvi.formats.agp import AGPLine

    line = AGPLine(line)
    if should_assert:
        with pytest.raises(AssertionError):
            line.validate()
    else:
        line.validate()
