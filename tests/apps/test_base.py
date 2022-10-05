#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import time

import pytest


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


def test_download():
    from jcvi.apps.base import cleanup, download
    from jcvi.apps.vecscreen import ECOLI_URL, UNIVEC_URL

    ret = download("http://www.google.com")
    assert ret == "index.html"
    cleanup(ret)

    ret = download(ECOLI_URL, filename="ecoli.fa.gz")
    assert ret == "ecoli.fa.gz"
    cleanup(ret)

    ret = download(UNIVEC_URL, filename="univec.fa.gz")
    assert ret == "univec.fa.gz"
    cleanup(ret)

    ret = download(UNIVEC_URL)
    assert ret == "UniVec_Core"
    cleanup(ret)


def test_ls_ftp():
    from jcvi.apps.base import ls_ftp

    url = "ftp://ftp.ensembl.org/pub/release-75/fasta/"
    valid_species = [x for x in ls_ftp(url) if "." not in x]
    assert "saccharomyces_cerevisiae" in valid_species
    assert "gorilla_gorilla" in valid_species
    assert len(valid_species) == 67


@pytest.mark.parametrize(
    "input_list,output_list",
    [
        ([], []),
        ([1], [1]),
        ([1, [2]], [1, 2]),
        ([1, [2, "33"]], [1, 2, "33"]),
        ([1, [2, "33"], ("45",)], [1, 2, "33", "45"]),
        ([[["aaa"], "bbb"], [[[["ccc"]]]]], ["aaa", "bbb", "ccc"]),
    ],
)
def test_flatten(input_list, output_list):
    from jcvi.apps.base import flatten

    assert flatten(input_list) == output_list


def test_cleanup():
    from jcvi.apps.base import cleanup
    from jcvi.formats.base import write_file

    write_file("a", "content_a", skipcheck=True)
    write_file("b", "content_b", skipcheck=True)
    write_file("c", "content_c", skipcheck=True)
    cleanup("a", ["b", "c"])
    assert not op.exists("a")
    assert not op.exists("b")
    assert not op.exists("c")


def test_need_update():
    from jcvi.apps.base import cleanup, need_update
    from jcvi.formats.base import write_file

    cleanup("a", "b", "c")
    assert need_update("does_not_exist.txt", "does_not_exist.txt")

    write_file("a", "content_a", skipcheck=True)
    assert not need_update("a", "a")

    time.sleep(0.1)
    write_file("b", "content_b", skipcheck=True)
    assert need_update("b", "a")
    assert not need_update("a", "b")

    time.sleep(0.1)
    write_file("c", "content_c", skipcheck=True)
    assert need_update("c", ["a", "b"])
    assert need_update(["c", "b"], "a")
    assert not need_update("a", ["b", "c"])
    assert need_update(["a", "b"], ["c", "d"])
    cleanup("a", "b", "c")
