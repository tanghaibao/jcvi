#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest

from jcvi.utils.validator import ValidationError
from jcvi.graphics.synteny import LayoutLine


@pytest.mark.parametrize(
    "row,delimiter",
    [
        ("0.5, 0.6,        0, center,    top,      ,     1,       grape Chr1", ","),
        (
            "0.5, 0.4,        0, leftalign, center, #fc8d62,     1, peach scaffold_1",
            ",",
        ),
        ("0.5, 0.6, 0, left, center, g", ","),
        ("0.25, 0.7, 45, center, center, m", ","),
    ],
)
def test_valid_layoutline(row, delimiter):
    assert LayoutLine(row, delimiter) is not None


@pytest.mark.parametrize(
    "row,delimiter,error",
    [
        (
            "0.5, 0.6,        0, center,    top,      ,     1,       grape Chr1",
            ".",
            ValueError,
        ),
        (
            "0.5, 0.4,        0, topalign, center, #fc8d62,     1, peach scaffold_1",
            ",",
            ValidationError,
        ),
        ("0.5, 0.6, 0, left, rightalign, g", ",", ValidationError),
        ("1.5, 0.7, 45, center, center, m", ",", ValidationError),
    ],
)
def test_invalid_layoutline(row, delimiter, error):
    with pytest.raises(error):
        _ = LayoutLine(row, delimiter)
