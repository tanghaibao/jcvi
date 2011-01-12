#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
find . -type l | %prog

Linux commands `touch` wouldn't modify the mtime for symlinks, this script can.
Use find to pipe in all the symlinks.
"""

import os
import os.path as op
import sys

from optparse import OptionParser


if __name__ == '__main__':

    p = OptionParser(__doc__)
    fp = sys.stdin

    for link_name in fp:
        link_name = link_name.strip()
        if not op.islink(link_name): continue

        # re-link the symlinks (similar to `ln -sf`)
        source = os.readlink(link_name)
        os.unlink(link_name)
        os.symlink(source, link_name)

