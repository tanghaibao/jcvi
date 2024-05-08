#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import pytest

from jcvi.apps.base import cleanup
from jcvi.formats.bed import merge
from jcvi.graphics.synteny import LayoutLine, main as synteny_main
from jcvi.utils.validator import ValidationError


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


def test_main():
    cwd = os.getcwd()
    os.chdir(op.join(op.dirname(__file__), "data"))
    cleanup("grape_peach.bed", "blocks.pdf")
    merge(["grape.bed", "peach.bed", "-o", "grape_peach.bed"])
    image_name = synteny_main(["blocks", "grape_peach.bed", "blocks.layout"])
    assert op.exists(image_name)
    os.chdir(cwd)
