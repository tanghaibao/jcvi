#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest

from jcvi.graphics.base import latex, markup, rc


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
        ("grape_grape vs peach_peach", r"grape\_grape vs peach\_peach"),
    ],
)
def test_latex(s, expected):
    assert latex(s) == expected, "Expect {}".format(expected)


def test_markup():
    rc("text", usetex=True)
    s = "Prupe_1G289800.1"
    assert markup(s) == r"Prupe\_1G289800.1"
    rc("text", usetex=False)
    assert markup(s) == s


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


def test_assign_array_missing_attribute():
    """Test that assign_array works when the attribute is not set on the object.

    This guards against AttributeError: 'LayoutLine' object has no attribute 'color'
    when an object without the attribute is present in the layout list.
    """
    from jcvi.graphics.base import AbstractLayout

    class SimpleItem:
        pass

    class MockLayout(AbstractLayout):
        def __init__(self):
            super().__init__("/dev/null")

    ml = MockLayout()
    item = SimpleItem()
    ml.append(item)

    # Should not raise AttributeError even if 'color' is not set on item
    ml.assign_array("color", ["red"])
    assert item.color == "red"


def test_assign_array_existing_attribute_preserved():
    """Test that assign_array does not overwrite an existing truthy attribute."""
    from jcvi.graphics.base import AbstractLayout

    class SimpleItem:
        def __init__(self, color):
            self.color = color

    class MockLayout(AbstractLayout):
        def __init__(self):
            super().__init__("/dev/null")

    ml = MockLayout()
    item_with_color = SimpleItem("blue")
    item_without_color = SimpleItem("")
    ml.append(item_with_color)
    ml.append(item_without_color)

    ml.assign_array("color", ["red", "green"])
    # Existing truthy color should be preserved
    assert item_with_color.color == "blue"
    # Empty/falsy color should be replaced
    assert item_without_color.color == "green"
