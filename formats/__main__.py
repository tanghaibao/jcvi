#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import os.path as op
from jcvi.apps.base import dmain


if __name__ == '__main__':
    cwd = op.dirname(__file__)
    dmain(cwd)
