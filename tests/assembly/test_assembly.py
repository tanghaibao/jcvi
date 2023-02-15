#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_assembly.py
# tests
#
# Created by Haibao Tang on 11/25/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#
from ..config import test_script, generate_tests


def pytest_generate_tests(metafunc):
    generate_tests(metafunc, "assembly")
