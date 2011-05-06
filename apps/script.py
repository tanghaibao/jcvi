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
    """
    %prog scriptname.py

    create a minimal boilerplate for a new script
    """
    p = OptionParser(main.__doc__)

    opts, args = p.parse_args()
    if len(args) != 1:
        sys.exit(p.print_help())

    script = args[0]

    if op.exists(script):
        logging.error("File `{0}` exists, remove it first".format(script))
    else:
        fw = open(script, "w")
        fw.write(template)
        fw.close()
        logging.debug("Template writes to `{0}`".format(script))


if __name__ == '__main__':
    main()
