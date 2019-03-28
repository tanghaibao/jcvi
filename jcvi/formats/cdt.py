#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog data.cdt data.nwk

Convert the result from Eisen's CLUSTER program: data.gtr and data.cdt into NEWICK format
"""

import sys
import csv
import logging

from collections import namedtuple
from itertools import groupby

from jcvi.apps.base import OptionParser
from jcvi.formats.base import BaseFile


GTRLine = namedtuple("GTRLine", "parent left_child right_child dist")


class CDT (BaseFile):

    def __init__(self, filename):
        super(CDT, self).__init__(filename)

        pf = filename.rsplit(".", 1)[0]
        self.gtrfile = pf + ".gtr"
        self.atrfile = pf + ".atr"
        self.get_names()

    def get_names(self):
        cdt_file = self.filename
        reader = csv.reader(file(cdt_file), delimiter="\t")

        gid = next(reader)
        assert gid[0] == "GID"
        aid = next(reader)
        if aid[0] == "AID":
            eweight = next(reader)
        else:
            eweight = aid
        assert eweight[0] == "EWEIGHT"

        self.gnames = [x[:2] for x in reader]
        self.anames = zip(aid, gid)[4:]

    def get_gtr_tree(self):

        from ete2 import Tree

        fp = open(self.gtrfile)
        reader = csv.reader(fp, delimiter="\t")
        nodes = {}
        gnames = dict(self.gnames)
        for g in map(GTRLine._make, reader):
            node = Tree()
            parent_name, parent_dist = g.parent, float(g.dist)
            for child in (g.left_child, g.right_child):
                if child in gnames:
                    node.add_child(name=gnames[child], dist=1 - parent_dist)
                else:
                    assert child in nodes, child
                    child_node, child_dist = nodes[child]
                    node.add_child(child_node, dist=child_dist - parent_dist)

            nodes[parent_name] = (node, parent_dist)

        self.gtr_tree = node

    def print_newick(self, nwk_file, gtr=True):

        self.gtr_tree.write(format=5, outfile=nwk_file)
        logging.debug("Newick tree written to `{0}`".format(nwk_file))

    def iter_partitions(self, cutoff=.3, gtr=True):
        from jcvi.utils.grouper import Grouper

        if gtr:
            names = self.gnames
            fp = open(self.gtrfile)
        else:
            names = self.anames
            fp = open(self.atrfile)

        reader = csv.reader(fp, delimiter="\t")
        grouper = Grouper()
        for g in map(GTRLine._make, reader):
            d = float(g.dist)
            if d < cutoff:
                continue

            grouper.join(g.parent, g.left_child, g.right_child)

        parents = {}
        for i, group in enumerate(grouper):
            for g in group:
                parents[g] = i

        partitions = [[parents.get(a, x), x] for a, x in names]
        for key, parts in groupby(partitions, key=lambda x: x[0]):
            yield list(x[1] for x in parts)


def main(args):

    cdt_file, nwk_file = args
    cdt = CDT(cdt_file)
    cdt.get_gtr_tree()
    cdt.print_newick(nwk_file)


if __name__ == '__main__':

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) != 2:
        sys.exit(not p.print_help())

    main(args)
