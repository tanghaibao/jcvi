#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Simple validator to make sure certain values match expectation.
"""

from typing import Collection, Union, TypeVar

ComparableType = Union[int, float]
T = TypeVar("T")


class ValidationError(Exception):
    pass


def validate_in_choices(value: T, choices: Collection[T], tag: str = "Value") -> bool:
    """
    Validate if certain value is among a collection.
    Args:
        value: value of interest
        choices (Collection): a collection (list, tuple, dict, set etc.) of values
        tag (str): the semantic meaning of value to be shown in error

    Returns:
        True if validation passes. Raises ValidationError if it fails
    """
    if value not in choices:
        raise ValidationError(f"{tag} must be one of {choices}, you have: {value}")
    return True


def validate_in_range(
    value: ComparableType,
    min_value: ComparableType,
    max_value: ComparableType,
    tag: str = "Value",
) -> bool:
    """
    Validate if certain value is numerically within range.

    Args:
        value: value of interest
        min_value: minimum expected value
        max_value: maximum expected value
        tag (str): the semantic meaning of value to be shown in error

    Returns:
        True if validation passes. Raises ValidationError if it fails.
    """
    if not min_value <= value <= max_value:
        raise ValidationError(
            f"{tag} must be between [{min_value}, {max_value}], you have: {value}"
        )
    return True
