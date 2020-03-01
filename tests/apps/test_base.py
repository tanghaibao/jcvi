#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os


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


def test_download():
    from jcvi.apps.base import download
    from jcvi.apps.vecscreen import ECOLI_URL, UNIVEC_URL

    ret = download("http://www.google.com")
    assert ret == "index.html"
    os.remove(ret)

    ret = download(ECOLI_URL, filename="ecoli.fa.gz")
    assert ret == "ecoli.fa.gz"
    os.remove(ret)

    ret = download(UNIVEC_URL, filename="univec.fa.gz")
    assert ret == "univec.fa.gz"
    os.remove(ret)

    ret = download(UNIVEC_URL)
    assert ret == "UniVec_Core"
    os.remove(ret)
