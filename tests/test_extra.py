#!/usr/bin/env python
# -*- coding: UTF-8 -*-


def test_apps_console():
    """ Test colorful ASCII text
    """
    from jcvi.apps.console import test
    test()


def test_formats_cblast():
    """ Test formats.cblast - BLAST parsing
    """
    from jcvi.formats.cblast import BlastLine

    b = BlastLine("Os09g11510	Os08g13650	92.31	39	3	0	2273	2311	3237	3199	0.001	54.0")
    assert b.query == 'Os09g11510'
    assert b.hitlen == 39
