#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Data structure to represent overlaps between pairwise sequences
"""

import sys

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('app', ''),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def app(args):
    """
    %prog app

    """
    p = OptionParser(app.__doc__)
    opts, args = p.parse_args(args)


if __name__ == '__main__':
    main()
