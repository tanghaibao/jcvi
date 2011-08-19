#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Connect to JCVI sybase account
"""

import os.path as op
import sys

from optparse import OptionParser

from jcvi.apps.base import mkdir, gethostname, getusername
from jcvi.apps.base import ActionDispatcher, sh, debug
debug()


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
        if "password" in row:
            password = _(row)

    hostname = gethostname()
    username = getusername()

    return hostname, username, password


def connect(db, sql):
    import Sybase

    db = Sybase.connect(hostname, username, password, db)
    c = db.cursor()
    c.execute(sql)
    return c.fetchall()


def main():

    actions = (
        ('libs', 'get list of lib_ids to to run by pull'),
        ('pull', 'pull the sequences from the TIGR database'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def libs(args):
    """
    %prog libs libfile

    Get list of lib_ids to be run by pull(). The SQL commands:

    select library.lib_id, library.name from library join bac on
        library.bac_id=bac.id where bac.lib_name="Medicago";
    select seq_name from sequence where seq_name like 'MBE%'
        and trash is null;
    """
    p = OptionParser(libs.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    libfile, = args

    db = "track"
    sqlcmd = "select library.lib_id, library.name, bac.gb# from library join bac on " + \
             "library.bac_id=bac.id where bac.lib_name='Medicago'"
    results = connect(db, sqlcmd)

    fw = open(libfile, "w")
    for lib_id, name, gb in results:
        name = name.translate(None, "\n")
        if not gb:
            gb = "None"

        print >> fw, "|".join((lib_id, name, gb))
    fw.close()


def pull(args):
    """
    %prog pull libfile

    Pull the sequences using the first column in the libfile.
    """
    p = OptionParser(pull.__doc__)
    p.add_option("--frag", default=False, action="store_true",
            help="The command to pull sequences from db [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    libfile, = args

    frag = opts.frag
    fp = open(libfile)
    hostname, username, password = get_profile()

    for row in fp:
        lib_id, name = row.split("|", 1)
        sqlfile = lib_id + ".sql"

        if not op.exists(sqlfile):
            fw = open(sqlfile, "w")
            print >> fw, "select seq_name from sequence where seq_name like" + \
                    " '{0}%' and trash is null".format(lib_id)
            fw.close()

        if frag:
            cmd = "pullfrag -D mtg2 -n {0}.sql -o {0} -q -S SYBPROD".format(lib_id)
            cmd += " -U {0} -P {1}".format(username, password)
        else:
            cmd = "pullseq -D mtg2 -n {0}.sql -o {0} -q".format(lib_id)
        sh(cmd)


if __name__ == '__main__':
    main()
