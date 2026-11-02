#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from subprocess import CalledProcessError

import pytest


@pytest.mark.parametrize("returncode,killed", [(137, True), (-9, True), (1, False)])
def test_log_lastal_error(returncode, killed, monkeypatch):
    import jcvi.apps.align as align

    messages = []
    monkeypatch.setattr(align.logger, "error", lambda msg, *a: messages.append(msg))
    align.log_lastal_error(CalledProcessError(returncode, "lastal", output=b"Killed"))
    assert any("out of memory" in x for x in messages) == killed
