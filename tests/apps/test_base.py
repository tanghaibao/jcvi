#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import mock
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
    from jcvi.apps.base import cleanup, mkdir
    from jcvi.formats.base import write_file

    write_file("a", "content_a", skipcheck=True)
    write_file("b", "content_b", skipcheck=True)
    write_file("c", "content_c", skipcheck=True)
    paths = ("a", "b", "c")
    for path in paths:
        assert op.exists(path)
    cleanup("a", ["b", "c"])
    for path in paths:
        assert not op.exists(path)

    # Test cleanup with a directory
    mkdir("adir")
    mkdir("bdir")
    write_file("bs", "content_bs", skipcheck=True)
    paths = ("adir", "bdir", "bs")
    for path in paths:
        assert op.exists(path)
    cleanup("adir", ["bdir", "bs"])
    for path in paths:
        assert not op.exists(path)


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


def test_set_image_options():
    from jcvi.apps.base import OptionParser

    with pytest.raises(Exception):
        p_fails = OptionParser(__doc__)
        p_fails.add_argument("--cov", default="jcvi", help="pytest coverage")
        p_fails.add_argument("--dpi", default=300, type=int, help="DPI")
        p_fails.set_image_options()

    # This should be fine
    p = OptionParser(__doc__)
    p.add_argument("--cov", default="jcvi", help="pytest coverage")
    p.add_argument("--seed", default=300, type=int, help="seed")
    p.add_argument("-s", default="dummy", action="store_true", help="dummy")

    p.set_image_options()

    # Try adding an argument group
    g = p.add_argument_group("group")
    g.add_argument("--group", default="jcvi", help="pytest coverage")


def test_getpath():
    from jcvi.apps.base import getpath

    with mock.patch("builtins.input", lambda _: "e\rvvvvvvvv = zzzzzzzz\n"):
        assert getpath("not-part-of-path", name="CLUSTALW2", warn="warn") is None

    with mock.patch("builtins.input", lambda _: "/bin"):
        test_cfg = "test.cfg"
        assert getpath("rm", name="rm", cfg=test_cfg, warn="warn") in (
            "/bin/rm",
            "/usr/bin/rm",
        )
