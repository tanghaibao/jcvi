#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Procedure to touch and copy softlinks
"""

import os
import os.path as op
import logging
import sys

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()

def main():
    
    actions = (
        ('touch', 'touch all the symlinks'),
        ('cp', 'cp all the symlinks to current folder'),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def touch(args):
    """
    find . -type l | %prog touch

    Linux commands `touch` wouldn't modify the mtime for symlinks, this script can.
    Use find to pipe in all the symlinks.
    """
    p = OptionParser(touch.__doc__)
    fp = sys.stdin

    for link_name in fp:
        link_name = link_name.strip()
        if not op.islink(link_name): continue
        if not op.exists(link_name): continue

        # re-link the symlinks (similar to `ln -sf`)
        source = os.readlink(link_name)
        os.unlink(link_name)
        os.symlink(source, link_name)


def cp(args):
    """
    find folder -type l | %prog cp
    """
    p = OptionParser(cp.__doc__)
    fp = sys.stdin

    for link_name in fp:
        link_name  = link_name.strip()
        if not op.islink(link_name): continue
        if not op.exists(link_name): continue

        source = os.readlink(link_name)
        link_dir = op.dirname(link_name)
        source = op.normpath(op.join(link_dir, source))
        source = op.abspath(source)
        link_name = op.basename(link_name)
        if not op.exists(link_name):
            os.symlink(source, link_name)
        logging.debug(" => ".join((source, link_name)))


if __name__ == '__main__':
    main()

