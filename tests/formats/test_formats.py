#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from ..config import test_script, generate_tests


def pytest_generate_tests(metafunc):
    generate_tests(metafunc, "formats")
