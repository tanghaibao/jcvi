#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest
from .config import test_script, generate_tests


def pytest_generate_tests(metafunc):
    generate_tests(metafunc, "compara")
