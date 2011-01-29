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

from jcvi.utils.cbook import human_size 
from jcvi.apps.base import ActionDispatcher, debug
debug()

def main():
    
    actions = (
        ('touch', 'touch all the symlinks'),
        ('cp', 'cp all the symlinks to current folder'),
        ('size', 'print the file sizes for the files pointed by symlinks'),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_abs_path(link_name):
    source = os.readlink(link_name)
    link_dir = op.dirname(link_name)
    source = op.normpath(op.join(link_dir, source))
    source = op.abspath(source)
    return source


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

        source = get_abs_path(link_name)

        # re-link the symlinks (similar to `ln -sf`)
        os.unlink(link_name)
        os.symlink(source, link_name)


def cp(args):
    """
    find folder -type l | %prog cp

    Copy all the softlinks to the current folder, using absolute paths
    """
    p = OptionParser(cp.__doc__)
    fp = sys.stdin

    for link_name in fp:
        link_name  = link_name.strip()
        if not op.islink(link_name): continue
        if not op.exists(link_name): continue

        source = get_abs_path(link_name)

        link_name = op.basename(link_name)
        if not op.exists(link_name):
            os.symlink(source, link_name)
        logging.debug(" => ".join((source, link_name)))


def size(args):
    """
    find folder -type l | %prog size 

    Get the size for all the paths that are pointed by the links
    """
    p = OptionParser(size.__doc__)
    fp = sys.stdin

    results = []
    for link_name in fp:
        link_name  = link_name.strip()
        if not op.islink(link_name): continue

        source = get_abs_path(link_name)

        link_name = op.basename(link_name)
        filesize = op.getsize(source)
        results.append((filesize, link_name))

    # sort by descending file size
    for filesize, link_name in sorted(results, reverse=True):
        filesize = human_size(filesize) 
        print >>sys.stderr, "%10s\t%s" % (filesize, link_name)


if __name__ == '__main__':
    main()

