#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import sys
import os.path as op
import logging

from optparse import OptionParser

from jcvi.formats.base import must_open
from jcvi.apps.base import ActionDispatcher, debug
debug()

default_template = """#!/usr/bin/env python
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

    if len(args) != 1:
        sys.exit(not p.print_help())


if __name__ == '__main__':
    main()
"""

graphic_template = """#!/usr/bin/env python
# -*- coding: UTF-8 -*-

\"\"\"
%prog datafile

Illustrate blablabla...
\"\"\"


import sys
import logging

from optparse import OptionParser

from jcvi.graphics.base import plt, set_image_options, savefig
from jcvi.apps.base import debug
debug()


def main():
    p = OptionParser(__doc__)
    opts, args, iopts = set_image_options(p)

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
"""


def main():
    """
    %prog scriptname.py

    create a minimal boilerplate for a new script
    """
    p = OptionParser(main.__doc__)
    p.add_option("-g", dest="graphic", default=False, action="store_true",
            help="create boilerplate for a graphic script")

    opts, args = p.parse_args()
    if len(args) != 1:
        sys.exit(p.print_help())

    script, = args
    template = graphic_template if opts.graphic else default_template

    fw = must_open(script, "w", checkexists=True)
    fw.write(template)
    fw.close()
    message = "template writes to `{0}`".format(script)
    if opts.graphic:
        message = "graphic " + message
    message = message[0].upper() + message[1:]
    logging.debug(message)


if __name__ == '__main__':
    main()
