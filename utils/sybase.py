#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Connect to JCVI sybase account
"""

import os.path as op
import sys
import logging

import re

from jcvi.apps.base import OptionParser

from jcvi.formats.base import must_open
from jcvi.apps.base import mkdir, gethostname, getusername
from jcvi.apps.base import ActionDispatcher, sh, debug
debug()


def db_defaults():
    """
    JCVI legacy Sybase database connection defaults
    """
    return "SYBPROD", "access", "access"


def get_profile(sqshrc="~/.sqshrc", hostname=None, username=None, password=None):
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

    dhost, duser, dpass = db_defaults()

    if not password:
        username, password = duser, dpass
    elif not username:
        username = getusername()

    if not hostname:
        hostname = dhost

    return hostname, username, password


def connect(dbname, hostname=None, username=None, password=None):
    import Sybase

    hostname, username, password = \
            get_profile(hostname=hostname, username=username, password=password)

    dbh = Sybase.connect(hostname, username, password, database=dbname)
    cur = dbh.cursor()
    return dbh, cur


def fetchall(cur, sql):
    cur.execute(sql)
    return cur.fetchall()


def execute(cur, sql):
    cur.execute(sql)


def commit(dbh):
    return dbh.commit()


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
    p = OptionParser(libs.__doc__)
    p.set_db_opts(dbname="track", credentials=None)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    libfile, = args

    sqlcmd = "select library.lib_id, library.name, bac.gb# from library join bac on " + \
             "library.bac_id=bac.id where bac.lib_name='Medicago'"
    cur = connect(opts.dbname)
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
    p = OptionParser(pull.__doc__)
    p.set_db_opts(dbname="mtg2", credentials=None)
    p.add_option("--frag", default=False, action="store_true",
            help="The command to pull sequences from db [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    libfile, = args

    dbname = opts.dbname
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
            cmd = "pullfrag -D {0} -n {1}.sql -o {1} -q -S {2}".format(dbname, lib_id, hostname)
            cmd += " -U {0} -P {1}".format(username, password)
        else:
            cmd = "pullseq -D {0} -n {1}.sql -o {1} -q".format(dbname, lib_id)
        sh(cmd)


to_commit_re = re.compile("|".join("^{0}".format(x) for x in ("update", "insert", "delete")), re.I)
def to_commit(query):
    """
    check if query needs to be committed (only if "update", "insert" or "delete")
    """
    if re.search(to_commit_re, query):
        return True
    return None


def query(args):
    """
    %prog query "SELECT feat_name FROM asm_feature WHERE feat_type = \\"{0}\\" AND end5 <= \\"{1}\\" AND end3 >= \\"{2}\\"" ::: datafile1 ....

    Script takes the data from tab-delimited datafile(s) and replaces the placeholders
    in the query which is then executed. Depending upon the type of query, results are
    either printed out (when running `select`) or not (when running `insert`, `update`
    or `delete`)

    If the query contains quotes around field values, then these need to be to be escaped with \\
    """
    p = OptionParser(query.__doc__)
    p.set_db_opts()
    p.add_option("--dryrun", default=False, action="store_true",
                 help="Don't commit to database. Just print queries [default: %default]")
    p.add_option("--fieldsep", default="\t",
                 help="Specify output field separator [default: '%default']")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

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
                atoms = row.strip().split("\t")
                mi = [int(x.strip("{}")) for x in m]
                assert max(mi) <= len(atoms), \
                        "Number of columns in `datafile`({0})".format(len(atoms)) + \
                        " != number of `placeholders`({0})".format(len(m))
                natoms = [atoms[x] for x in mi]
                nqry = qry
                for idx, (match, atom) in enumerate(zip(m, natoms)):
                    nqry = nqry.replace(match, atom)
                queries.add(nqry)
    else:
        queries.add(qry)

    if not opts.dryrun:
        fw = must_open(opts.outfile, "w")
        dbh, cur = connect(opts.dbname, hostname=opts.hostname, username=opts.username,\
                password=opts.password)
    cflag = None
    for qry in queries:
        if opts.dryrun:
            print qry
        else:
            if to_commit(qry):
                execute(cur, qry)
                cflag = True
            else:
                results = fetchall(cur, qry)
                for result in results:
                    print >> fw, fieldsep.join([str(x) for x in result])
    if not opts.dryrun and cflag:
        commit(dbh)


if __name__ == '__main__':
    main()
