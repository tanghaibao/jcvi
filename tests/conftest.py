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
    Copy a test data directory into the isolated working dir.

    Copy rather than symlink: outputs whose names collide with inputs
    (`grape_peach.bed`, `calibrate.json`) would write through to the checkout.
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

    Many entry points write relative to the cwd; autouse keeps new tests
    isolated without each one opting in.
    """
    monkeypatch.chdir(tmp_path)
