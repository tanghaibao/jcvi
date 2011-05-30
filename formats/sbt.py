#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
GenBank submission template. In nested blocks below.
"""

import csv
import sys
import string
import logging

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()

NameTemplate = """        {{
          name name {{
            last "{last}",
            first "{first}",
            initials "{initials}",
            suffix "{suffix}"
          }}
        }}"""

class SubmissionTemplate (object):
    def __init__(self):
        raise NotImplementedError


def main():

    actions = (
        ('names', 'convert a list of names to sbt blocks'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def names(args):
    """
    %prog names namelist templatefile

    Generate name blocks by accepting namelist file. namelist file is
    tab-delimited that contains >=4 columns of data. 
    
    Three columns are mandatory. First name, middle initiial and last name.
    First row is table header.

    For the extra columns, the first column will go in the `$N0` field
    in the template file, second to the `$N1` field, etc.
    """
    p = OptionParser(names.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    namelist, templatefile = args
    reader = csv.reader(open(namelist), delimiter="\t")
    header = reader.next()
    ncols = len(header)
    assert ncols > 3
    nextras = ncols - 3

    blocks = []
    bools = []
    for row in reader:
        first, middle, last = row[:3]
        extras = row[3:]
        bools.append([(x.upper() == 'Y') for x in extras])
        middle = middle.strip()
        if middle != "":
            middle = middle.rstrip('.') + '.'
        initials = "{0}.{1}".format(first[0], middle)
        suffix = ""
        nameblock = NameTemplate.format(last=last, first=first,
                initials=initials, suffix=suffix)
        blocks.append(nameblock)

    selected_idx = zip(*bools)
    out = [] * nextras
    for i, sbools in enumerate(selected_idx):
        selected = []
        for b, ss in zip(blocks, sbools):
            if ss:
                selected.append(b)
        bigblock = ",\n".join(selected)
        out.append(bigblock)
        logging.debug("List N{0} contains a total of {1} names.".format(i,
            len(selected)))

    template = open(templatefile).read()
    template = string.Template(template)
    outmapping = dict(("N{0}".format(i), x) for (i, x) in enumerate(out))
    t = template.substitute(outmapping)
    print t


if __name__ == '__main__':
    main()
