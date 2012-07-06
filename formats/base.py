#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import math
import sys
import logging

from itertools import groupby, islice, cycle, izip
from optparse import OptionParser

from Bio import SeqIO
from jcvi.apps.base import ActionDispatcher, sh, debug, need_update, \
        mkdir, popen, set_outfile
debug()


class BaseFile (object):

    def __init__(self, filename):

        self.filename = filename
        if filename:
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
        ncols = max(keypos, valuepos) + 1
        thiscols = 0
        for lineno, row in enumerate(fp):
            row = row.rstrip()
            atoms = row.split(delimiter)
            thiscols = len(atoms)
            if thiscols < ncols:
                msg = "Must contain >= {0} columns.  Aborted.\n".format(ncols)
                msg += "  --> Line {0}: {1}".format(lineno + 1, row)
                logging.error(msg)
                sys.exit(1)

            key = atoms[keypos]
            value = atoms[valuepos] if (valuepos is not None) else atoms
            self[key] = value

        assert thiscols, "File empty"
        self.ncols = thiscols
        logging.debug("Imported {0} records from `{1}`.".\
                    format(len(self), filename))


class SetFile (BaseFile, set):

    def __init__(self, filename, column=-1, delimiter=None):
        super(SetFile, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            keys = [x.strip() for x in row.split(delimiter)]
            if column >= 0:
                keys = [keys[column]]
            self.update(keys)


class FileShredder (object):
    """
    Same as rm -f *
    """
    def __init__(self, filelist):

        cmd = "rm -f {0}".format(" ".join(filelist))
        sh(cmd)


class FileMerger (object):
    """
    Same as cat * > filename
    """
    def __init__(self, filelist, outfile):

        self.filelist = filelist
        self.outfile = outfile

    def merge(self, checkexists=False):
        outfile = self.outfile
        if checkexists and not need_update(self.filelist, outfile):
            logging.debug("File `{0}` exists. Merge skipped.".format(outfile))
            return

        files = " ".join(self.filelist)
        sh("cat {0}".format(files), outfile=outfile)


class FileSplitter (object):

    def __init__(self, filename, outputdir=None, mode="cycle"):
        self.filename = filename
        self.outputdir = outputdir
        self.mode = mode

        self.format = format = self._guess_format(filename)
        logging.debug("format is %s" % format)

        if format in ("fasta", "fastq"):
            self.klass = "seqio"
        else:
            self.klass = "txt"

        mkdir(outputdir)

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

        if ext in ("fasta", "fa", "fna", "cds", "pep", "faa"):
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
        root, ext = op.splitext(op.basename(filename))

        names = []
        for i in xrange(N):
            name = "%s_%02d%s" % (root, i, ext)
            names.append(name)

        return names

    def split(self, N, force=False):
        """
        There are two modes of splitting the records
        - batch: splitting is sequentially to records/N chunks
        - cycle: placing each record in the splitted files and cycles

        use `cycle` if the len of the record is not evenly distributed
        """
        mode = self.mode
        assert mode in ("batch", "cycle")
        logging.debug("set split mode=%s" % mode)

        self.names = self.__class__.get_names(self.filename, N)
        if self.outputdir:
            self.names = [op.join(self.outputdir, x) for x in self.names]

        if not need_update(self.filename, self.names) and not force:
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


def check_exists(filename):
    """
    Avoid overwriting some files accidentally.
    """
    if op.exists(filename):
        logging.error("`{0}` found, overwrite (Y/N)?".format(filename))
        overwrite = (raw_input() == 'Y')
    else:
        overwrite = True

    return overwrite


def must_open(filename, mode="r", checkexists=False, skipcheck=False):
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

    elif filename == "tmp" and mode == "w":
        from tempfile import NamedTemporaryFile
        fp = NamedTemporaryFile(delete=False)

    elif filename.endswith(".gz"):
        if 'r' in mode:
            cmd = "zcat {0}".format(filename)
            fp = popen(cmd, debug=False)
        elif 'w' in mode:
            import gzip
            fp = gzip.open(filename, mode)

    elif filename.endswith(".bz2"):
        if 'r' in mode:
            cmd = "bzcat {0}".format(filename)
            fp = popen(cmd, debug=False)
        elif 'w' in mode:
            import bz2
            fp = bz2.BZ2File(filename, mode)

    else:
        if checkexists:
            assert mode == "w"
            overwrite = (not op.exists(filename)) if skipcheck \
                        else check_exists(filename)
            if overwrite:
                fp = open(filename, "w")
            else:
                logging.debug("File `{0}` already exists. Skipped."\
                        .format(filename))
                return None
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
        ('reorder', 'reorder columns in tab-delimited files'),
        ('flatten', 'convert a list of IDs into one per line'),
        ('setop', 'set operations on files'),
        ('join', 'join tabular files based on common column'),
        ('truncate', 'remove lines from end of file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def truncate(args):
    """
    %prog truncate linecount filename

    Remove linecount lines from the end of the file in-place. Borrowed from:
    <http://superuser.com/questions/127786/how-to-remove-the-last-2-lines-of-a-very-large-file>
    """
    p = OptionParser(truncate.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    number, filename = args
    number = int(number)
    count = 0

    f = open(filename, "r+b")
    f.seek(0, os.SEEK_END)
    end = f.tell()
    while f.tell() > 0:
        f.seek(-1, os.SEEK_CUR)
        char = f.read(1)
        if char == '\n':
            count += 1
        if count == number + 1:
            f.truncate()
            print >> sys.stderr, "Removed {0} lines from end of file".format(number)
            return number

        f.seek(-1, os.SEEK_CUR)

    if count < number + 1:
        print >> sys.stderr, "No change: requested removal would leave empty file"
        return -1


def flatten(args):
    """
    %prog flatten filename > ids

    Convert a list of IDs (say, multiple IDs per line) and move them into one
    per line.
    """
    p = OptionParser(flatten.__doc__)
    p.add_option("--sep", default=",",
                 help="Separater for the tabfile [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    tabfile, = args
    fp = open(tabfile)
    for row in fp:
        print row.strip().replace(opts.sep, "\n")


def reorder(args):
    """
    %prog reorder tabfile 1,2,4,3 > newtabfile

    Reorder columns in tab-delimited files. The above syntax will print out a
    new file with col-1,2,4,3 from the old file.
    """
    import csv

    p = OptionParser(reorder.__doc__)
    p.add_option("--sep", default="\t",
                 help="Separater for the tabfile [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    tabfile, order = args
    sep = opts.sep
    order = [int(x) - 1 for x in order.split(",")]
    reader = csv.reader(must_open(tabfile), delimiter=sep)
    writer = csv.writer(sys.stdout, delimiter=sep)
    for row in reader:
        newrow = [row[x] for x in order]
        writer.writerow(newrow)


def split(args):
    """
    %prog split file outdir

    split file into records
    """
    p = OptionParser(split.__doc__)
    p.add_option("-n", dest="N", type="int", default=1,
            help="split into N chunks [default: %default]")
    p.add_option("--all", default=False, action="store_true",
            help="split all records [default: %default]")
    p.add_option("--cycle", default=False, action="store_true",
            help="splitted records in Round Robin fashion [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    mode = "cycle" if opts.cycle else "batch"
    filename, outdir = args
    fs = FileSplitter(filename, outputdir=outdir, mode=mode)

    if opts.all:
        logging.debug("option -all override -n")
        N = fs.num_records
    else:
        N = opts.N

    logging.debug("split file into %d chunks" % N)
    fs.split(N)

    return fs


def join(args):
    """
    %prog join file1.txt file2.txt ..

    Join tabular files based on common column. --column specifies the column
    index to pivot on. Use comma to separate multiple values if the pivot column
    is different in each file. Maintain the order in the first file.
    """
    from jcvi.utils.iter import flatten

    p = OptionParser(join.__doc__)
    p.add_option("--column", default="0",
                 help="0-based column id, multiple values allowed [default: %default]")
    p.add_option("--noheader", default=False, action="store_true",
                 help="Do not print header [default: %default]")
    set_outfile(p)

    opts, args = p.parse_args(args)
    nargs = len(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    c = opts.column
    if "," in c:
        cc = [int(x) for x in c.split(",")]
    else:
        cc = [int(c)] * nargs

    assert len(cc) == nargs, "Column index number != File number"

    # Maintain the first file line order, and combine other files into it
    pivotfile = args[0]
    files = [DictFile(f, keypos=c, valuepos=None, delimiter="\t") \
                        for f, c in zip(args, cc)]
    otherfiles = files[1:]
    header = "\t".join(flatten([op.basename(x.filename)] * x.ncols \
                        for x in files))

    fp = open(pivotfile)
    fw = must_open(opts.outfile, "w")
    if not opts.noheader:
        print >> fw, header

    for row in fp:
        row = row.rstrip()
        atoms = row.split("\t")
        newrow = atoms
        key = atoms[cc[0]]
        for d in otherfiles:
            drow = d.get(key, ["na"] * d.ncols)
            newrow += drow
        print >> fw, "\t".join(newrow)


def setop(args):
    """
    %prog setop "fileA & fileB" > newfile

    Perform set operations, except on files. The files (fileA and fileB) contain
    list of ids. The operator is one of the four:

    |: union (elements found in either file)
    &: intersection (elements found in both)
    -: difference (elements in fileA but not in fileB)
    ^: symmetric difference (elementes found in either set but not both)

    Please quote the argument to avoid shell interpreting | and &.
    """
    p = OptionParser(setop.__doc__)
    p.add_option("--column", default=0, type="int",
                 help="The column to extract, 0-based, -1 to disable [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    statement, = args
    fa, op, fb = statement.split()
    assert op in ('|', '&', '-', '^')

    column = opts.column
    fa = SetFile(fa, column=column)
    fb = SetFile(fb, column=column)

    if op == '|':
        t = fa | fb
    elif op == '&':
        t = fa & fb
    elif op == '-':
        t = fa - fb
    elif op == '^':
        t = fa ^ fb

    for x in sorted(t):
        print x


if __name__ == '__main__':
    main()
