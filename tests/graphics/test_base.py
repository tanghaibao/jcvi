#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest


@pytest.mark.parametrize(
    "s,expected",
    [
        ("a short name", "a short name"),
        (
            "a really really long name for you a really really long name for you",
            "a really...e for you",
        ),
        ("These colors look lovely together", "These co... together"),
    ],
)
def test_shorten(s, expected):
    from jcvi.graphics.base import shorten

    assert shorten(s) == expected, "Expect {}".format(expected)


@pytest.mark.parametrize(
    "s,expected",
    [
        ("grape_grape vs peach_peach", "grape\_grape vs peach\_peach"),
    ],
)
def test_latex(s, expected):
    from jcvi.graphics.base import latex

    assert latex(s) == expected, "Expect {}".format(expected)


@pytest.mark.parametrize(
    "figname,format,expected",
    [
        ("my_test.pdf", "pdf", "my_test.pdf"),
        ("my_test.png", "pdf", "my_test.png"),
        ("sss", "pdf", "sss.pdf"),
        ("aaa.bbb.svg", "pdf", "aaa.bbb.svg"),
        ("xxx.yyy", "svg", "xxx.yyy.svg"),
    ],
)
def test_update_figname(figname, format, expected):
    from jcvi.graphics.base import update_figname

    assert update_figname(figname, format) == expected, "Expect {}".format(expected)
