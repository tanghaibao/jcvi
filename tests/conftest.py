#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Shared pytest fixtures for the jcvi test suite.
"""

import os
import os.path as op
import shutil

import pytest


@pytest.fixture
def copy_data(tmp_path):
    """
    Copy the contents of a test data directory into the isolated working dir.

    Returns a callable taking the data directory to copy from.

    Tests must *copy* rather than symlink: several entry points write outputs
    whose names collide with their inputs (`grape_peach.bed`, `calibrate.json`),
    and writing through a symlink modifies the link target, i.e. the checked-in
    data directory. Copying keeps those writes inside the temporary directory.
    """

    def _copy(data_dir: str) -> str:
        for fname in os.listdir(data_dir):
            src, dst = op.join(data_dir, fname), op.join(str(tmp_path), fname)
            if op.isdir(src):
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)
        return str(tmp_path)

    return _copy


@pytest.fixture(autouse=True)
def isolate_cwd(tmp_path, monkeypatch):
    """
    Run every test in its own temporary directory.

    Many jcvi entry points write their outputs relative to the current working
    directory. Without this, running the suite from the repo root leaves stray
    files behind in the checkout and in `tests/*/data` directories. Making the
    isolation automatic means new tests get it for free, rather than each one
    having to remember to `monkeypatch.chdir(tmp_path)`.

    Tests that need the real data directories should reference them via
    `op.join(op.dirname(__file__), "data")`, which is unaffected by the chdir.
    """
    monkeypatch.chdir(tmp_path)
