import os
import os.path as op
import math
import sys
import logging

from itertools import groupby, islice, cycle, izip
from optparse import OptionParser

from Bio import SeqIO
from jcvi.apps.base import ActionDispatcher, sh, debug
debug()


class BaseFile (object):

    def __init__(self, filename):

        self.filename = filename
        logging.debug("Load file `{0}`".format(filename))


class LineFile (BaseFile, list):
    """
    Generic file parser for line-based files
    """
    def __init__(self, filename):

        super(LineFile, self).__init__(filename)


class DictFile (BaseFile, dict):
    """
    Generic file parser for multi-column files, keyed by a particular index.
    """
    def __init__(self, filename, keypos=0, valuepos=1, delimiter=None):

        super(DictFile, self).__init__(filename)

        fp = must_open(filename)
        for lineno, row in enumerate(fp):
            row = row.rstrip()
            atoms = row.split(delimiter)
            if len(atoms) < 2:
                msg = "Line must contain >= 2 columns. Aborted.\n"
                msg += "  --> Line {0}: {1}".format(lineno + 1, row)
                logging.error(msg)
                sys.exit(1)

            key, value = atoms[keypos], atoms[valuepos]
            self[key] = value


class FileMerger (object):
    """
    Same as cat * > filename
    """
    def __init__(self, filelist, outfile):

        self.filelist = filelist
        self.outfile = outfile

    def merge(self):
        files = " ".join(self.filelist)
        sh("cat {0}".format(files), outfile=self.outfile)


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

        if self.klass == "seqio":
            handle = SeqIO.parse(open(filename), self.format)
        else:
            handle = open(filename)
        return handle

    @property
    def num_records(self):
        handle = self._open(self.filename)
        return sum(1 for x in handle)

    def _guess_format(self, filename):
        root, ext = op.splitext(filename)
        ext = ext.strip(".")

        if ext in ("fasta", "fa", "fna", "cds", "pep"):
            format = "fasta"
        elif ext in ("fastq",):
            format = "fastq"
        else:
            format = "txt"
        return format

    def _batch_iterator(self, N=1):
        """Returns N lists of records.

        This can be used on any iterator, for example to batch up
        SeqRecord objects from Bio.SeqIO.parse(...), or to batch
        Alignment objects from Bio.AlignIO.parse(...), or simply
        lines from a file handle.

        This is a generator function, and it returns lists of the
        entries from the supplied iterator.  Each list will have
        batch_size entries, although the final list may be shorter.
        """
        batch_size = math.ceil(self.num_records / float(N))
        handle = self._open(self.filename)
        while True:
            batch = list(islice(handle, batch_size))
            if not batch:
                break
            yield batch

    @classmethod
    def get_names(cls, filename, N):
        root, ext = op.splitext(filename)

        names = []
        for i in xrange(N):
            name = "%s_%02d%s" % (root, i, ext)
            names.append(name)

        return names

    def split(self, N, force=False, mode="cycle"):
        """
        There are two modes of splitting the records
        - batch: splitting is sequentially to records/N chunks
        - cycle: placing each record in the splitted files and cycles

        use `cycle` if the len of the record is not evenly distributed
        """
        assert mode in ("batch", "cycle")
        logging.debug("set split mode=%s" % mode)

        self.names = self.__class__.get_names(self.filename, N)
        if self.outputdir:
            self.names = [op.join(self.outputdir, x) for x in self.names]

        if op.exists(self.names[0]) and not force:
            logging.error("file %s already existed, skip file splitting" % \
                    self.names[0])
            return

        filehandles = [open(x, "w") for x in self.names]

        if mode == "batch":
            for batch, fw in zip(self._batch_iterator(N), filehandles):

                if self.klass == "seqio":
                    count = SeqIO.write(batch, fw, self.format)
                else:
                    for line in batch:
                        fw.write(line)
                    count = len(batch)

                logging.debug("write %d records to %s" % (count, fw.name))

        elif mode == "cycle":
            handle = self._open(self.filename)
            for record, fw in izip(handle, cycle(filehandles)):

                if self.klass == "seqio":
                    SeqIO.write(record, fw, self.format)
                else:
                    fw.write(record)

        for fw in filehandles:
            fw.close()


def must_open(filename, mode="r"):
    """
    Accepts filename and returns filehandle.

    Checks on multiple files, stdin/stdout/stderr, .gz or .bz2 file.
    """
    if isinstance(filename, list):
        assert "r" in mode

        import fileinput
        return fileinput.input(filename)

    if filename in ("-", "stdin"):
        assert "r" in mode
        fp = sys.stdin

    elif filename == "stdout":
        assert "w" in mode
        fp = sys.stdout

    elif filename == "stderr":
        assert "w" in mode
        fp = sys.stderr

    elif filename.endswith(".gz"):
        import gzip
        fp = gzip.open(filename, mode)

    elif filename.endswith(".bz2"):
        import bz2
        fp = bz2.BZ2File(filename, mode)

    else:
        fp = open(filename, mode)

    return fp


def read_until(handle, start):
    # read each line until a certain start, then puts the start tag back
    while 1:
        pos = handle.tell()
        line = handle.readline()
        if not line:
            break
        if line.startswith(start):
            handle.seek(pos)
            return
    #raise EOFError, "%s tag cannot be found"


def read_block(handle, signal):
    """
    Useful for reading block-like file formats, for example FASTA or OBO file,
    such file usually startswith some signal, and in-between the signals are a
    record
    """
    signal_len = len(signal)
    it = (x[1] for x in groupby(handle,
        key=lambda row: row.strip()[:signal_len] == signal))
    for header in it:
        header = header.next().strip()
        if header[:signal_len] != signal:
            continue
        seq = list(s.strip() for s in it.next())
        yield header, seq


def main():

    actions = (
        ('split', 'split large file into N chunks'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def split(args):
    """
    %prog split file outdir

    split file into records
    """
    p = OptionParser(split.__doc__)
    p.add_option("-n", dest="N", type="int", default=1,
            help="split into N chunks")
    p.add_option("--all", default=False, action="store_true",
            help="split all records")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    filename, outdir = args
    fs = FileSplitter(filename, outputdir=outdir)

    if opts.all:
        logging.debug("option -all override -n")
        N = fs.num_records
    else:
        N = opts.N

    logging.debug("split file into %d chunks" % N)
    fs.split(N)


if __name__ == '__main__':
    main()
