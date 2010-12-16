
"""
Classes to handle the .bed filed$ and .raw file
"""

import sys

from base import LineFile


class BedLine(object):
    # the Bed format supports more columns. we only need
    # the first 4, but keep the information in 'stuff'.
    __slots__ = ("seqid", "start", "end", "accn", "stuff")

    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.seqid = args[0]
        self.start = int(args[1])
        self.end = int(args[2])
        self.accn = args[3]
        self.stuff = args[4:] if len(args) > 4 else None

    def __str__(self):
        s = "\t".join(map(str, [getattr(self, attr) \
                    for attr in BedLine.__slots__[:-1]]))
        if self.stuff:
            s += "\t" + "\t".join(self.stuff)
        return s

    def __getitem__(self, key):
        return getattr(self, key)


class Bed(LineFile):

    def __init__(self, filename=None, key=None):
        super(Bed, self).__init__(filename)

        # the sorting key provides some flexibility in ordering the features
        # for example, user might not like the lexico-order of seqid
        self.key = key or (lambda x: (x.seqid, x.start, x.accn))

        if not filename: return

        for line in open(filename):
            if line[0] == "#": continue
            if line.startswith('track'): continue
            self.append(BedLine(line))

        self.sort(key=self.key)

    def print_to_file(self, fw=sys.stdout):
        for bedline in self:
            print >>fw, bedline

    @property
    def seqids(self):
        return sorted(set(b.seqid for b in self))

    @property
    def order(self):
        # get the gene order given a Bed object
        return dict((f.accn, (i, f)) for (i, f) in enumerate(self))

    @property
    def simple_bed(self):
        return [(b.seqid, i) for (i, b) in enumerate(self)]


