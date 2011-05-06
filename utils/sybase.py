#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Connect to JCVI sybase account
"""

import os.path as op


def get_profile(sqshrc="~/.sqshrc"):
    """
    get database, username, password from .sqshrc file e.g.
    \set username="user"
    """
    _ = lambda x: x.split("=")[-1].translate(None, "\"'").strip()
    sqshrc = op.expanduser(sqshrc)
    for row in open(sqshrc):
        row = row.strip()
        if not row.startswith("\set"):
            continue
        if "username" in row:
            username = _(row)
        elif "password" in row:
            password = _(row)
        elif "hostname" in row:
            hostname = _(row)

    return hostname, username, password

hostname, username, password = get_profile()


def connect(db, sql):
    import Sybase

    db = Sybase.connect(hostname, username, password, db)
    c = db.cursor()
    c.execute(sql)
    return c.fetchall()


if __name__ == '__main__':
    print get_profile()
    sql = "select bac_name from bac where bac_name like 'mth2%'"
    print connect("track", sql)
