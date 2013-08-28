#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Connect to JCVI sybase account
"""

import os.path as op
import sys
import logging

from jcvi.apps.base import MOptionParser

from jcvi.formats.base import must_open
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
        if not row.startswith("\\set") or "prompt" in row:
            continue
        if "password" in row:
            password = _(row)
        if "hostname" in row:
            hostname = _(row)
        if "username" in row:
            username = _(row)

    if not username:
        username = getusername()
    if not hostname:
        hostname = gethostname()

    return hostname, username, password


def connect(db):
    import Sybase

    hostname, username, password = get_profile()
    db = Sybase.connect(hostname, username, password, database=db)
    c = db.cursor()
    return c


def fetchall(cur, sql):
    cur.execute(sql)
    return cur.fetchall()


def runsql(cur, sql):
    return cur.execute(sql)


def main():

    actions = (
        ('libs', 'get list of lib_ids to to run by pull'),
        ('pull', 'pull the sequences from the TIGR database'),
        ('query', 'run query using input from datafile'),
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
    p = MOptionParser(libs.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    libfile, = args

    db = "track"
    sqlcmd = "select library.lib_id, library.name, bac.gb# from library join bac on " + \
             "library.bac_id=bac.id where bac.lib_name='Medicago'"
    cur = connect(db)
    results = fetchall(cur, sqlcmd)

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
    p = MOptionParser(pull.__doc__)
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


def query(args):
    """
    %prog query "SELECT feat_name FROM asm_feature
        WHERE feat_type = \\"{0}\\"
        AND end5 <= \\"{1}\\"
        AND end3 >= \\"{2}\\"" ::: datafile(s)

    Script takes the data from tab-delimited `datafile` and replaces the placeholders
    in the query which is then executed. Depending upon the type of query, results are
    either printed out (when running select) or not (when running insert, update, delete)

    If the query contains quotes around field values, then these need to be to be escaped with \\
    """
    import re

    valid_qtypes = ["select", "insert", "update", "delete"]

    p = MOptionParser(query.__doc__)
    p.add_option("--db", default="mta4",
                 help="Specify name of database to query [default: %default]")
    p.add_option("--qtype", default="select", choices=valid_qtypes,
                 help="Specify type of query being run [default: %default]")
    p.add_option("--fieldsep", default="\t",
                 help="Specify output field separator [default: '%default']")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    dbname = opts.db
    qtype = opts.qtype
    fieldsep = opts.fieldsep

    sep = ":::"
    files = None
    if sep in args:
        sepidx = args.index(sep)
        files = args[sepidx + 1:]
        args = args[:sepidx]
        if not files:
            files = [""]

    qry = " ".join(args)
    queries = set()
    if files:
        m = re.findall(r"\{\d+\}", qry)
        if not m:
            logging.error("Error: query `{0}` does not contain placeholders".format(qry))
            sys.exit()
        for datafile in files:
            datafile = datafile.strip()
            fp = must_open(datafile)
            for row in fp:
                atoms = row.split()
                assert len(atoms) == len(m), \
                        "Number of columns in `datafile`({0})".format(len(atoms)) + \
                        " != number of `placeholders`({0})".format(len(m))
                mi = [int(x.strip("{}")) for x in m]
                natoms = [x for (mi,x) in sorted(zip(mi,atoms))]
                for idx, (match, atom) in enumerate(zip(m, natoms)):
                    nqry = qry.replace(match, atom)
                    queries.add(nqry)
    else:
        queries.add(qry)

    fw = must_open(opts.outfile, "w")
    cur = connect(dbname)
    for qry in queries:
        if qtype == "select":
            results = fetchall(cur, qry)
            for result in results:
                print >> fw, fieldsep.join([str(x) for x in result])
        else:
            status = runsql(cur, qry)
            if status is not None:
                logging.debug("Error in query {0}".format(qry))


if __name__ == '__main__':
    main()
