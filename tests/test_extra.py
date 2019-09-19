#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest


class TestAppsConsole(unittest.TestCase):
    """ Test colorful ASCII text
    """

    def test_apps_console(self):
        from jcvi.apps.console import test

        test()


class TestFormatsCblast(unittest.TestCase):
    """ Test formats.cblast - BLAST parsing
    """

    def setUp(self):
        from jcvi.formats.cblast import BlastLine

        # Test parser
        self.b = BlastLine(
            "Os09g11510	Os08g13650	92.31	39	3	0	2273	2311	3237	3199	0.001	54.0"
        )

    def test_formats_getters(self):
        assert self.b.query == "Os09g11510"
        assert self.b.hitlen == 39

    def test_formats_setters(self):
        # Test setters
        self.b.query = "At5g54690"
        assert str(self.b).startswith("At5g54690")

    def test_formats_swap(self):
        query, subject = self.b.query, self.b.subject
        self.b.query, self.b.subject = subject, query
        assert str(self.b).startswith("Os08g13650")
