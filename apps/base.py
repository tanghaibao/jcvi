"""
basic support for running library as script
"""

import os
import os.path as op
import sys
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

    if opts.debug: debug()


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


def sh(cmd, grid=False, outfile=None):
    """
    simple wrapper for system calls
    """
    if grid:
        from jcvi.apps.grid import GridProcess
        pr = GridProcess(cmd, outfile=outfile)
        pr.start(path=None)
    else:
        cmd += " > {0}".format(outfile)
        logging.debug(cmd)
        call(cmd, shell=True)


def is_current_file(a, b):
    """
    Check if the file a is newer than file b
    """
    if not (op.exists(a) and op.exists(b)): return False
    am = os.stat(a).st_mtime
    bm = os.stat(b).st_mtime
    return am > bm


def get_today():
    """
    Returns the date in 2010-07-14 format
    """
    from datetime import date
    return str(date.today())


def debug():
    """
    turn on the debugging
    """
    import logging
    
    from jcvi.apps.console import ColoredText

    format = "%(asctime)s [%(module)s::%(levelname)s] %(message)s"
    format = str(ColoredText(format, "yellow"))
    logging.basicConfig(level=logging.DEBUG,
            format=format,
            datefmt="%H:%M:%S",
            )

