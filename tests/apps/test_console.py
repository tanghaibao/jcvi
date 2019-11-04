#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest


class TestConsole(unittest.TestCase):
    """ Test colorful ASCII text
    """

    def test_smoke(self):
        from jcvi.apps.console import test

        test()
