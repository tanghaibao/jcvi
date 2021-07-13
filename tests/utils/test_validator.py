#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import pytest
from jcvi.utils.validator import ValidationError, validate_in_choices, validate_in_range


@pytest.mark.parametrize("value,min_value,max_value", [(0.5, 0, 1), (250, 200, 300)])
def test_valid_in_range(value, min_value, max_value):
    assert validate_in_range(value, min_value, max_value)


@pytest.mark.parametrize("value,min_value,max_value", [(-0.5, 0, 1), (1000, 200, 300)])
def test_invalid_in_range(value, min_value, max_value):
    with pytest.raises(ValidationError):
        validate_in_range(value, min_value, max_value)


@pytest.mark.parametrize("value,choices", [(7, range(10)), ("aa", ["a", "aa", "bc"])])
def test_valid_in_choices(value, choices):
    assert validate_in_choices(value, choices)


@pytest.mark.parametrize("value,choices", [(10, range(10)), ("aaa", ["a", "aa", "bc"])])
def test_invalid_in_choices(value, choices):
    with pytest.raises(ValidationError):
        validate_in_choices(value, choices)
