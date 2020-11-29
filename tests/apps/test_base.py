#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op


def test_sample_N():
    from jcvi.apps.base import sample_N

    # No random seed, we test the length of the output
    a = sample_N(list(range(10)), 2)
    assert len(a) == 2
    assert set(a) < set(range(10))

    a = sample_N(list(range(10)), 10)
    assert len(a) == 10
    assert set(range(2)) < set(a)

    a = sample_N(list(range(5)), 5)
    assert len(a) == 5
    assert set(a) == set(range(5))

    # With random seed, we test for exact matches
    a = sample_N([1, 2, 3], 2, seed=666)
    assert a == [2, 3]

    a = sample_N([1, 2, 3], 3, seed=666)
    assert a == [2, 3, 1]

    a = sample_N([1, 2, 3], 4, seed=666)
    assert a == [2, 3, 1, 2]


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


def test_ls_ftp():
    from jcvi.apps.base import ls_ftp

    url = "ftp://ftp.ensembl.org/pub/release-75/fasta/"
    valid_species = [x for x in ls_ftp(url) if "." not in x]
    assert "saccharomyces_cerevisiae" in valid_species
    assert "gorilla_gorilla" in valid_species
    assert len(valid_species) == 67
