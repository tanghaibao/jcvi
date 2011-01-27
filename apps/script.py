#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import sys
import os.path as op
import logging

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()

template = """#!/usr/bin/env python
# -*- coding: UTF-8 -*-

\"\"\"

\"\"\"

import sys

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('app', ''),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def app(args):
    \"\"\"
    %prog app

    \"\"\"
    p = OptionParser(app.__doc__)
    opts, args = p.parse_args(args)


if __name__ == '__main__':
    main()
"""

def main():

    actions = (
        ('make', 'make boilerplate script'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def make(args):
    """
    %prog make scriptname.py

    create a minimal boilerplate for a new script
    """
    p = OptionParser(make.__doc__)
    p.add_option("-f", dest="force", default=False, action="store_true",
            help="force overwrite [default: %default]")    

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())
    
    script = args[0]

    if op.exists(script) and not opts.force:
        logging.debug("file `%s` exists, use -f to force overwrite" % script)
    else:
        fw = open(script, "w")
        fw.write(template)
        fw.close()
        logging.debug("template writes to file `%s`" % script)


if __name__ == '__main__':
    main()
