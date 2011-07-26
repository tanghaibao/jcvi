"""
basic support for running library as script
"""

import os
import os.path as op
import sys
import shutil
import logging

from subprocess import call
from optparse import OptionParser


class ActionDispatcher (object):

    def __init__(self, actions):

        self.actions = actions
        self.valid_actions, self.action_helps = zip(*actions)

    def print_help(self):
        help = "Available actions:\n"
        for action, action_help in self.actions:
            help += "    `%s`: %s\n" % (action, action_help)

        print >>sys.stderr, help
        sys.exit(1)

    def dispatch(self, globals):
        if len(sys.argv) == 1:
            self.print_help()

        action = sys.argv[1]

        if not action in self.valid_actions:
            print >>sys.stderr, "%s not a valid action" % action
            self.print_help()

        globals[action](sys.argv[2:])


def set_debug(instance, args):
    """
    Add --debug options for command line programs
    """
    assert isinstance(instance, OptionParser)

    instance.add_option("--debug", dest="debug",
            default=False, action="store_true",
            help="set to debug level")

    opts, args = instance.parse_args(args)

    if opts.debug:
        debug()


def set_grid(instance):
    """
    Add --grid options for command line programs
    """
    assert isinstance(instance, OptionParser)

    instance.add_option("--grid", dest="grid",
            default=False, action="store_true",
            help="run on the grid [default: %default]")


def set_params(instance):
    """
    Add --params options for given command line programs
    """
    assert isinstance(instance, OptionParser)

    instance.add_option("--params", dest="extra", default="",
            help="extra parameters to run")


def sh(cmd, grid=False, infile=None, outfile=None, errfile=None,
        background=False):
    """
    simple wrapper for system calls
    """
    if grid:
        from jcvi.apps.grid import GridProcess
        pr = GridProcess(cmd, infile=infile, outfile=outfile, errfile=errfile)
        pr.start(path=None)
        return 0  # A fake retcode
    else:
        if infile:
            cmd += " < {0} ".format(infile)
        if outfile:
            cmd += " > {0} ".format(outfile)
        if errfile:
            cmd += " 2> {0} ".format(errfile)
        if background:
            cmd += " & "

        logging.debug(cmd)
        return call(cmd, shell=True)


def popen(cmd):
    """
    Capture the cmd stdout output to a file handle.
    """
    from subprocess import Popen, PIPE
    logging.debug(cmd)
    proc = Popen(cmd, bufsize=1, stdout=PIPE, shell=True)
    return proc.stdout


def mkdir(dirname, overwrite=False):
    """
    Wraps around os.mkdir(), but checks for existence first.
    """
    if op.isdir(dirname):
        if overwrite:
            shutil.rmtree(dirname)
            os.mkdir(dirname)
            logging.debug("Overwrite folder `{0}`.".format(dirname))
        else:
            return False  # Nothing is changed
    else:
        os.mkdir(dirname)
        logging.debug("`{0}` not found. Creating new.".format(dirname))

    return True


def is_newer_file(a, b):
    """
    Check if the file a is newer than file b
    """
    if not (op.exists(a) and op.exists(b)):
        return False
    am = os.stat(a).st_mtime
    bm = os.stat(b).st_mtime
    return am > bm


def get_today():
    """
    Returns the date in 2010-07-14 format
    """
    from datetime import date
    return str(date.today())


def download(url, filename=None):
    from urlparse import urlsplit

    scheme, netloc, path, query, fragment = urlsplit(url)
    filename = filename or op.basename(path)

    if op.exists(filename):
        msg = "File `{0}` exists. Download skipped.".format(filename)
        logging.error(msg)
    else:
        cmd = "wget {0} -O {1}".format(url, filename)
        sh(cmd)

    return filename


def debug():
    """
    Turn on the debugging
    """
    import logging

    from jcvi.apps.console import ColoredText

    format = "%(asctime)s [%(module)s::%(levelname)s] %(message)s"
    format = str(ColoredText(format, "yellow"))
    logging.basicConfig(level=logging.DEBUG,
            format=format,
            datefmt="%H:%M:%S",
            )
