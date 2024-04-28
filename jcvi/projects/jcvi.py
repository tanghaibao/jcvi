#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Functions in this script produce figures in the JCVI manuscript.
"""

from ..apps.base import ActionDispatcher, OptionParser, logger


def genomebuild(args):
    """
    %prog genomebuild

    Plot genome build composite figure.
    """
    p = OptionParser(genomebuild.__doc__)
    _, args, iopts = p.set_image_options(args, figsize="12x9")


def main():

    actions = (
        ("genomebuild", "Plot genome build composite figure"),
        ("diversity", "Plot diversity composite figure"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
