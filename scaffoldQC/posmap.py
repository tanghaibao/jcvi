#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
POSMAP (POSitional MAPping) files are part of the  Celera Assembler output. 

Specs:
http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=POSMAP
"""

import sys
import csv

from collections import namedtuple
from optparse import OptionParser

from jcvi.formats.base import LineFile
from jcvi.apps.base import ActionDispatcher, debug
debug()

class PosmapFrags (object):
    pass


PosmapMatesLine = namedtuple("PosmapMatesLine", 
        "firstReadID secondReadID mateStatus")

class PosmapMates (object):

    def __init__(self, filename):
        fp = csv.reader(open(filename), delimiter='\t')
        for row in fp:
            b = PosmapMatesLine._make(row)
            print b


class Posmap (LineFile):

    # dispatch based on filename
    mapping = {
            "frags": PosmapFrags,
            "mates": PosmapMates,
            }

    def __init__(self, filename):
        super(Posmap, self).__init__(filename)


    def parse(self):
        # this is essentially to dispatch data
        filename = self.filename
        suffix = filename.rsplit(".", 1)[-1]
        assert suffix in self.mapping, \
                "`{0}` unknown format".format(filename)

        klass = self.mapping[suffix]
        return klass(filename)


def main():

    actions = (
        ('parse', 'parse posmap file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def parse(args):
    """
    %prog parse posmap

    """
    p = OptionParser(parse.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    posmapfile = args[0]
    p = Posmap(posmapfile)

    p.parse()


if __name__ == '__main__':
    main()
