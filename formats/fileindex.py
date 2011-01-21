"""
FileIndex
=========

index flat files. see: `blogpost <http://hackmap.blogspot.com/2010/04/fileindex.html>`_
example::

    >>> fi = FileIndex(f, FastQEntry)
    >>> print fi
    FileIndex(filename=`t.fastq`)
    >>> print ','.join(fi.keys()[:4])
    @SNPSTER3:1:1:4:1119#0/2,@SNPSTER3:1:1:4:1153#0/2,@SNPSTER3:1:1:4:1311#0/2,@SNPSTER3:1:1:4:1932#0/2
    >>> print fi["@SNPSTER3:1:1:4:1952#0/2"].id
    @SNPSTER3:1:1:4:1952#0/2
    >>> print fi[2].id
    @SNPSTER3:1:1:4:1311#0/2

adapted from brentp's fileindex and adapted to local convention
<https://github.com/brentp/bio-playground/tree/master/fileindex>
use bsddb instead of tcdb

This is meant as a generic solution to indexing flatfiles, if indexing sequence
files, use biopython's index instead
"""

import os
import sys
import os.path as op
import logging
import bsddb

from jcvi.apps.base import debug
from jcvi.formats.base import read_until
debug()

class FileIndex(object):
    ext = ".fidx"

    def __init__(self, filename, call_class, key=lambda x: x.id, 
            mode='c', allow_multiple=False):
        self.filename = filename
        self.allow_multiple = allow_multiple
        self.fh = open(self.filename)
        self.call_class = call_class
        self.key = key
        self.idxfile = filename + FileIndex.ext
        # only create and read are supported
        assert mode in ('c', 'r')

        if op.exists(self.idxfile) and mode=='c': 
            self.clear()

        self.db = bsddb.btopen(self.idxfile, mode)
        if mode=='c': self.create()

    def __getitem__(self, key):
        # supports indexing of both integer and string
        if isinstance(key, int):
            key = self.keys()[key]
        pos = self.db[key]
        self.fh.seek(long(pos))
        return self.call_class(self.fh)

    def __repr__(self):
        return "FileIndex(filename=`%s`)" % self.filename
    
    def create(self):
        fh = self.fh
        fh.seek(0)
        pos = fh.tell()
        while True:
            key = self.key(self.call_class(fh))
            if not key: break
            self.db[key] = str(pos) 
            #print "%s => %s" % (key, pos)
            # fh has been moved forward by get_next.
            pos = fh.tell()

    def close(self):
        self.fh.close()
        self.db.close()

    def clear(self):
        logging.debug("drop idex file `%s`" % self.idxfile)
        os.remove(self.idxfile)

    def keys(self):
        return self.db.keys()


class FastaEntry (object):
    def __init__(self, fh):
        line = fh.readline().strip()
        self.id = line.split()[0] if line else None
        read_until(fh, ">")


if __name__ == "__main__":

    f = sys.argv[1]
    fi = FileIndex(f, FastaEntry)
    print fi
    print ','.join(fi.keys()[:4])
    print fi[2].id
    fi.close()

