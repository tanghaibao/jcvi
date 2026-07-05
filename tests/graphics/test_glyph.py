#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from jcvi.graphics.glyph import demo


def test_demo(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    demo([])
