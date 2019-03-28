#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Procedure to touch and copy softlinks
"""
from __future__ import print_function

import os
import os.path as op
import logging
import sys

from jcvi.apps.base import OptionParser, ActionDispatcher, get_abs_path


def main():

    actions = (
        ('touch', 'touch all the symlinks'),
        ('cp', 'cp all the symlinks to current folder'),
        ('clean', 'removes all the symlinks in current folder'),
        ('size', 'print the file sizes for the files pointed by symlinks'),
        ('link', 'link source to target based on a tabular file'),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def lnsf(source, target, log=False):
    # re-link the symlinks (similar to `ln -sf`)
    if op.lexists(target):
        os.unlink(target)
    os.symlink(source, target)
    if log:
        logging.debug("{0} => {1}".format(source, target))


def link(args):
    """
    %prog link metafile

    Link source to target based on a tabular file.
    """
    from jcvi.apps.base import mkdir

    p = OptionParser(link.__doc__)
    p.add_option("--dir",
                 help="Place links in a subdirectory [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    meta, = args
    d = opts.dir
    if d:
        mkdir(d)

    fp = open(meta)
    cwd = op.dirname(get_abs_path(meta))
    for row in fp:
        source, target = row.split()
        source = op.join(cwd, source)
        if d:
            target = op.join(d, target)
        lnsf(source, target, log=True)


def touch(args):
    """
    find . -type l | %prog touch

    Linux commands `touch` wouldn't modify mtime for links, this script can.
    Use find to pipe in all the symlinks.
    """
    p = OptionParser(touch.__doc__)
    opts, args = p.parse_args(args)
    fp = sys.stdin

    for link_name in fp:
        link_name = link_name.strip()
        if not op.islink(link_name):
            continue
        if not op.exists(link_name):
            continue

        source = get_abs_path(link_name)
        lnsf(source, link_name)


def clean(args):
    """
    %prog clean

    Removes all symlinks from current folder
    """
    p = OptionParser(clean.__doc__)
    opts, args = p.parse_args(args)

    for link_name in os.listdir(os.getcwd()):
        if not op.islink(link_name):
            continue
        logging.debug("remove symlink `{0}`".format(link_name))
        os.unlink(link_name)


def cp(args):
    """
    find folder -type l | %prog cp

    Copy all the softlinks to the current folder, using absolute paths
    """
    p = OptionParser(cp.__doc__)
    fp = sys.stdin

    for link_name in fp:
        link_name = link_name.strip()
        if not op.exists(link_name):
            continue

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
    from jcvi.utils.cbook import human_size

    p = OptionParser(size.__doc__)
    fp = sys.stdin

    results = []
    for link_name in fp:
        link_name = link_name.strip()
        if not op.islink(link_name):
            continue

        source = get_abs_path(link_name)

        link_name = op.basename(link_name)
        filesize = op.getsize(source)
        results.append((filesize, link_name))

    # sort by descending file size
    for filesize, link_name in sorted(results, reverse=True):
        filesize = human_size(filesize, a_kilobyte_is_1024_bytes=True)
        print("%10s\t%s" % (filesize, link_name), file=sys.stderr)


if __name__ == '__main__':
    main()
