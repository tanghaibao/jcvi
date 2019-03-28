#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import sys
import logging

from jcvi.formats.base import write_file
from jcvi.apps.base import OptionParser


default_template = """
\"\"\"

\"\"\"

import sys

{}
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('app', ''),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

{}

if __name__ == '__main__':
    main()
"""

default_imports = ""

graphic_imports = """
from jcvi.graphics.base import plt, savefig, normalize_axes"""

default_app = """
def app(args):
    \"\"\"
    %prog app datafile

    \"\"\"
    p = OptionParser(app.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())
"""

graphic_app = """
def app(args):
    \"\"\"
    %prog app datafile

    \"\"\"
    p = OptionParser(app.__doc__)
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)
"""


def main():
    """
    %prog scriptname.py

    create a minimal boilerplate for a new script
    """
    p = OptionParser(main.__doc__)
    p.add_option("-g", "--graphic", default=False, action="store_true",
            help="Create boilerplate for a graphic script")

    opts, args = p.parse_args()
    if len(args) != 1:
        sys.exit(not p.print_help())

    script, = args
    imports = graphic_imports if opts.graphic else default_imports
    app = graphic_app if opts.graphic else default_app
    template = default_template.format(imports, app)
    write_file(script, template)

    message = "template writes to `{0}`".format(script)
    if opts.graphic:
        message = "graphic " + message
    message = message[0].upper() + message[1:]
    logging.debug(message)


if __name__ == '__main__':
    main()
