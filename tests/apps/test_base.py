#!/usr/bin/env python
# -*- coding: UTF-8 -*-


def test_sample_N():
    from jcvi.apps.base import sample_N

    a = sample_N(list(range(10)), 2)
    assert len(a) == 2
    assert set(a) < set(range(10))

    a = sample_N(list(range(10)), 10)
    assert len(a) == 10
    assert set(range(2)) < set(a)

    a = sample_N(list(range(5)), 5)
    assert len(a) == 5
    assert set(a) == set(range(5))
