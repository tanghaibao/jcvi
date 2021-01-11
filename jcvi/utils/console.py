#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# console.py
# utils
#
# Created by Haibao Tang on 01/09/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#

"""
We create a singleton console instance at the module level or as an attribute
of your top-level object.
"""

from rich.console import Console

console = Console()
printf = console.print
