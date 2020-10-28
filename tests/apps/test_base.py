#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op


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


def test_remove_if_exists():
    from jcvi.apps.base import remove_if_exists

    filename = "test_remove_if_exists.txt"
    remove_if_exists(filename)  # nothing happens

    with open(filename, "w") as fw:
        print("0", file=fw)

    assert op.exists(filename)
    remove_if_exists(filename)
    assert not op.exists(filename)


def test_download():
    from jcvi.apps.base import download, remove_if_exists
    from jcvi.apps.vecscreen import ECOLI_URL, UNIVEC_URL

    ret = download("http://www.google.com")
    assert ret == "index.html"
    remove_if_exists(ret)

    ret = download(ECOLI_URL, filename="ecoli.fa.gz")
    assert ret == "ecoli.fa.gz"
    remove_if_exists(ret)

    ret = download(UNIVEC_URL, filename="univec.fa.gz")
    assert ret == "univec.fa.gz"
    remove_if_exists(ret)

    ret = download(UNIVEC_URL)
    assert ret == "UniVec_Core"
    remove_if_exists(ret)
