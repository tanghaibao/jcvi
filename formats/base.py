import os
import os.path as op
import math
import sys
import logging
from optparse import OptionParser

from Bio import SeqIO
from jcvi.apps.base import ActionDispatcher
from jcvi.utils.cbook import sh

class BaseFile (object):
    
    def __init__(self, filename):

        self.filename = filename
        logging.debug("Load file %s" % filename)


class LineFile (BaseFile, list):
    """
    Generic file parser for line-based files
    """
    def __init__(self, filename):

        super(LineFile, self).__init__(filename)


class FileMerger (object):
    """
    same as cat * > filename
    """
    def __init__(self, filelist, outfile):
        
        self.filelist = filelist
        self.outfile = outfile

    def merge(self):
        files = " ".join(self.filelist)
        sh("cat %s > %s" % (files, self.outfile))


class FileSplitter (object):

    def __init__(self, filename, outputdir=None):
        self.filename = filename
        self.outputdir = outputdir

        self.format = format = self._guess_format(filename)
        logging.debug("format is %s" % format)

        if format in ("fasta", "fastq"):
            self.klass = "seqio"
        else:
            self.klass = "txt"

        if not op.isdir(outputdir):
            os.mkdir(outputdir)

    def _open(self, filename):

        if self.klass=="seqio":
            handle = SeqIO.parse(open(filename), self.format)
        else:
            handle = open(filename)
        return handle 

    def _num_records(self):
        handle = self._open(self.filename)
        return sum(1 for x in handle)

    def _guess_format(self, filename):
        root, ext = op.splitext(filename)
        ext = ext.strip(".")

        if ext in ("fasta", "fa", "fna"):
            format = "fasta"
        elif ext in ("fastq",):
            format = "fastq"
        else:
            format = "txt"
        return format

    def _batch_iterator(self, N=1) :
        """Returns N lists of records.
     
        This can be used on any iterator, for example to batch up
        SeqRecord objects from Bio.SeqIO.parse(...), or to batch
        Alignment objects from Bio.AlignIO.parse(...), or simply
        lines from a file handle.
     
        This is a generator function, and it returns lists of the
        entries from the supplied iterator.  Each list will have
        batch_size entries, although the final list may be shorter.
        """

        entry = True #Make sure we loop once
        batch_size = math.ceil(self._num_records() / float(N))
        handle = self._open(self.filename)
        while entry :
            batch = []
            while len(batch) < batch_size :
                try :
                    entry = handle.next()
                except StopIteration :
                    entry = None
                if entry is None: # end of file
                    break
                batch.append(entry)
            if batch:
                yield batch

    @classmethod
    def get_names(cls, filename, N):
        root, ext = op.splitext(filename)

        names = []
        for i in xrange(N):
            name = "%s_%02d%s" % (root, i, ext)
            names.append(name)

        return names

    def split(self, N, force=False):
        
        assert 1 < N <100, "number of pieces must be 1 < N < 100"
        self.names = self.__class__.get_names(self.filename, N)

        for batch, filename in zip(self._batch_iterator(N), self.names):
            if self.outputdir:
                filename = op.join(self.outputdir, filename)

            if op.exists(filename) and not force:
                logging.error("file %s already existed, skip file splitting" % \
                        filename)
                return

            fw = open(filename, "w")

            if self.klass=="seqio":
                count = SeqIO.write(batch, fw, self.format)
            else:
                for line in batch: fw.write(line)
                count = len(batch)

            fw.close()
            logging.debug("wrote %d records to %s" % (count, filename))


def main():
    
    actions = (
        ('split', 'split large file into N chunks'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def split(args):
    """
    %prog split file N 
    
    split file into N chunks
    """
    p = OptionParser(split.__doc__)

    logging.basicConfig(level=logging.DEBUG)
    opts, args = p.parse_args(args)

    try:
        filename, N = args
        N = int(N)
    except Exception, e:
        logging.error(str(e))
        sys.exit(p.print_help())

    fs = FileSplitter(filename)
    fs.split(N)


if __name__ == '__main__':
    main()
