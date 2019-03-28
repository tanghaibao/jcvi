#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Machine learning algorithms.
"""
from __future__ import print_function

import sys

from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('libsvm', 'convert csv file to LIBSVM format'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def libsvm(args):
    """
    %prog libsvm csvfile prefix.ids

    Convert csv file to LIBSVM format. `prefix.ids` contains the prefix mapping.
    Ga -1
    Gr 1

    So the feature in the first column of csvfile get scanned with the prefix
    and mapped to different classes. Formatting spec:

    http://svmlight.joachims.org/
    """
    from jcvi.formats.base import DictFile

    p = OptionParser(libsvm.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    csvfile, prefixids = args
    d = DictFile(prefixids)
    fp = open(csvfile)
    next(fp)
    for row in fp:
        atoms = row.split()
        klass = atoms[0]
        kp = klass.split("_")[0]
        klass = d.get(kp, "0")
        feats = ["{0}:{1}".format(i + 1, x) for i, x in enumerate(atoms[1:])]
        print(" ".join([klass] + feats))


if __name__ == '__main__':
    main()
