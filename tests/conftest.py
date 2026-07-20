#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Shared pytest fixtures for the jcvi test suite.
"""

import pytest


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
