#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
import os
import os.path as op
import math
import sys
import logging

from collections import OrderedDict
from itertools import groupby, islice, cycle

from Bio import SeqIO
from jcvi.apps.base import (
    OptionParser,
    ActionDispatcher,
    sh,
    debug,
    need_update,
    mkdir,
    popen,
)

debug()


FastaExt = ("fasta", "fa", "fna", "cds", "pep", "faa", "fsa", "seq", "nt", "aa")
FastqExt = ("fastq", "fq")


class BaseFile(object):
    def __init__(self, filename):

        self.filename = filename
        if filename:
            logging.debug("Load file `{0}`".format(filename))


class LineFile(BaseFile, list):
    """
    Generic file parser for line-based files
    """

    def __init__(self, filename, comment=None, load=False):

        super(LineFile, self).__init__(filename)

        if load:
            fp = must_open(filename)
            self.lines = [l.strip() for l in fp if l[0] != comment]
            logging.debug(
                "Load {0} lines from `{1}`.".format(len(self.lines), filename)
            )


class DictFile(BaseFile, OrderedDict):
    """
    Generic file parser for multi-column files, keyed by a particular index.
    """

    def __init__(
        self,
        filename,
        keypos=0,
        valuepos=1,
        delimiter=None,
        strict=True,
        keycast=None,
        cast=None,
    ):

        BaseFile.__init__(self, filename)
        OrderedDict.__init__(self)
        self.keypos = keypos

        fp = must_open(filename)
        ncols = (max(keypos, valuepos) if valuepos else keypos) + 1
        thiscols = 0
        for lineno, row in enumerate(fp):
            row = row.rstrip()
            atoms = row.split(delimiter)
            atoms = [x.strip() for x in atoms]
            thiscols = len(atoms)
            if thiscols < ncols:
                action = "Aborted" if strict else "Skipped"

                msg = "Must contain >= {0} columns.  {1}.\n".format(ncols, action)
                msg += "  --> Line {0}: {1}".format(lineno + 1, row)
                logging.error(msg)
                if strict:
                    sys.exit(1)
                else:
                    continue

            key = atoms[keypos]
            value = atoms[valuepos] if (valuepos is not None) else atoms
            if keycast:
                key = keycast(key)
            if cast:
                value = cast(value)
            self[key] = value

        assert thiscols, "File empty"
        self.ncols = thiscols
        logging.debug("Imported {0} records from `{1}`.".format(len(self), filename))

    @classmethod
    def num_columns(cls, filename, delimiter=None):
        """Return the column number of the csv file.

        Args:
            filename (str): Path to the file.
            delimiter (str, optional): Separator of the csv file. Defaults to None.

        Returns:
            int: Column number.
        """
        fp = must_open(filename)
        return max(len(row.split(delimiter)) for row in fp)


class SetFile(BaseFile, set):
    def __init__(self, filename, column=-1, delimiter=None):
        super(SetFile, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            if not row.strip():
                continue
            keys = [x.strip() for x in row.split(delimiter)]
            if column >= 0:
                keys = [keys[column]]
            self.update(keys)


class FileShredder(object):
    """
    Same as rm -f *
    """

    def __init__(self, filelist, verbose=True):

        filelist = [x for x in filelist if x and op.exists(x)]
        cmd = "rm -rf {0}".format(" ".join(filelist))
        sh(cmd, log=verbose)


class FileMerger(object):
    """
    Same as cat * > filename
    """

    def __init__(self, filelist, outfile):

        self.filelist = filelist
        self.outfile = outfile
        self.ingz = filelist[0].endswith(".gz")
        self.outgz = outfile.endswith(".gz")

    def merge(self, checkexists=False):
        outfile = self.outfile
        if checkexists and not need_update(self.filelist, outfile):
            logging.debug("File `{0}` exists. Merge skipped.".format(outfile))
            return

        files = " ".join(self.filelist)
        ingz, outgz = self.ingz, self.outgz
        if ingz and outgz:  # can merge gz files directly
            cmd = "cat {0} > {1}".format(files, outfile)
            sh(cmd)
        else:
            cmd = "zcat" if self.ingz else "cat"
            cmd += " " + files
            sh(cmd, outfile=outfile)

        return outfile


class FileSplitter(object):
    def __init__(self, filename, outputdir=None, format="fasta", mode="cycle"):
        self.filename = filename
        self.outputdir = outputdir
        self.mode = mode

        format = format or self._guess_format(filename)
        logging.debug("format is %s" % format)

        if format in ("fasta", "fastq"):
            self.klass = "seqio"
        elif format == "clust":
            self.klass = "clust"
        else:
            self.klass = "txt"

        self.format = format
        mkdir(outputdir)

    def _open(self, filename):

        if self.klass == "seqio":
            handle = SeqIO.parse(open(filename), self.format)
        elif self.klass == "clust":
            from jcvi.apps.uclust import ClustFile

            handle = iter(ClustFile(filename))
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

        if ext in FastaExt:
            format = "fasta"
        elif ext in FastqExt:
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
        pad0 = len(str(int(N - 1)))
        for i in range(N):
            name = "{0}_{1:0{2}d}{3}".format(root, i, pad0, ext)
            names.append(name)

        return names

    def write(self, fw, batch):
        if self.klass == "seqio":
            SeqIO.write(batch, fw, self.format)
        elif self.klass == "clust":
            for b in batch:
                print(b, file=fw)
        else:
            for line in batch:
                fw.write(line)
        return len(batch)

    def split(self, N, force=False):
        """
        There are two modes of splitting the records
        - batch: splitting is sequentially to records/N chunks
        - cycle: placing each record in the splitted files and cycles

        use `cycle` if the len of the record is not evenly distributed
        """
        mode = self.mode
        assert mode in ("batch", "cycle", "optimal")
        logging.debug("set split mode=%s" % mode)

        self.names = self.__class__.get_names(self.filename, N)
        if self.outputdir:
            self.names = [op.join(self.outputdir, x) for x in self.names]

        if not need_update(self.filename, self.names) and not force:
            logging.error(
                "file %s already existed, skip file splitting" % self.names[0]
            )
            return

        filehandles = [open(x, "w") for x in self.names]

        if mode == "batch":
            for batch, fw in zip(self._batch_iterator(N), filehandles):
                count = self.write(fw, batch)
                logging.debug("write %d records to %s" % (count, fw.name))

        elif mode == "cycle":
            handle = self._open(self.filename)
            for record, fw in zip(handle, cycle(filehandles)):
                count = self.write(fw, [record])

        elif mode == "optimal":
            """
            This mode is based on Longest Processing Time (LPT) algorithm:

            A simple, often-used algorithm is the LPT algorithm (Longest
            Processing Time) which sorts the jobs by its processing time and
            then assigns them to the machine with the earliest end time so far.
            This algorithm achieves an upper bound of 4/3 - 1/(3m) OPT.

            Citation: <http://en.wikipedia.org/wiki/Multiprocessor_scheduling>
            """
            endtime = [0] * N
            handle = self._open(self.filename)
            for record in handle:
                mt, mi = min((x, i) for (i, x) in enumerate(endtime))
                fw = filehandles[mi]
                count = self.write(fw, [record])
                endtime[mi] += len(record)

        for fw in filehandles:
            fw.close()


def longest_unique_prefix(query, targets, remove_self=True):
    """
    Find the longest unique prefix for filename, when compared against a list of
    filenames. Useful to simplify file names in a pool of files. See usage in
    formats.fasta.pool().
    """
    query = op.basename(query)
    targets = [op.basename(x) for x in targets]
    prefix_lengths = [len(op.commonprefix([query, name])) for name in targets]
    if remove_self and len(query) in prefix_lengths:
        prefix_lengths.remove(len(query))
    longest_length = max(prefix_lengths)
    return query[: longest_length + 1]


def check_exists(filename, oappend=False):
    """
    Avoid overwriting some files accidentally.
    """
    from jcvi.utils.console import console

    if op.exists(filename):
        if oappend:
            return oappend
        overwrite = (
            console.input("`{}` found, overwrite (Y/n)?".format(filename)) == "Y"
        )
    else:
        overwrite = True

    return overwrite


def timestamp():
    from datetime import datetime as dt

    return "{0}{1:02d}{2:02d}".format(dt.now().year, dt.now().month, dt.now().day)


def must_open(filename, mode="r", checkexists=False, skipcheck=False, oappend=False):
    """
    Accepts filename and returns filehandle.

    Checks on multiple files, stdin/stdout/stderr, .gz or .bz2 file.
    """
    if isinstance(filename, list):
        assert "r" in mode

        if filename[0].endswith((".gz", ".bz2")):
            filename = " ".join(filename)  # allow opening multiple gz/bz2 files
        else:
            import fileinput

            return fileinput.input(filename)

    if filename.startswith("s3://"):
        from jcvi.utils.aws import pull_from_s3

        filename = pull_from_s3(filename)

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

        fp = NamedTemporaryFile(mode=mode, delete=False)

    elif filename.endswith(".gz"):
        import gzip

        if "r" in mode:
            fp = gzip.open(filename, mode + "t")
        elif "w" in mode:
            fp = gzip.open(filename, mode)

    elif filename.endswith(".bz2"):
        if "r" in mode:
            cmd = "bzcat {0}".format(filename)
            fp = popen(cmd, debug=False)
        elif "w" in mode:
            import bz2

            fp = bz2.BZ2File(filename, mode)

    else:
        if checkexists:
            assert mode == "w"
            overwrite = (
                (not op.exists(filename))
                if skipcheck
                else check_exists(filename, oappend)
            )
            if overwrite:
                if oappend:
                    fp = open(filename, "a")
                else:
                    fp = open(filename, "w")
            else:
                logging.debug("File `{0}` already exists. Skipped.".format(filename))
                return None
        else:
            fp = open(filename, mode)

    return fp


bash_shebang = "#!/bin/bash"
python_shebang = """#!/usr/bin/env python
# -*- coding: UTF-8 -*-"""


def write_file(filename, contents, meta=None, skipcheck=False, append=False, tee=False):
    if not meta:
        suffix = filename.rsplit(".", 1)[-1]
        if suffix == "sh":
            meta = "run script"
        elif suffix == "py":
            meta = "python script"
        else:
            meta = "file"

    meta_choices = ("file", "run script", "python script")
    assert meta in meta_choices, "meta must be one of {0}".format(
        "|".join(meta_choices)
    )

    contents = contents.strip()
    shebang = "\n"
    if "script" in meta:
        if not append:
            if meta == "run script":
                shebang = bash_shebang
            elif meta == "python script":
                shebang = python_shebang
        contents = "\n\n".join((shebang, contents))

    fw = must_open(filename, "w", checkexists=True, skipcheck=skipcheck, oappend=append)
    if fw:
        print(contents, file=fw)
        fw.close()
    if tee:
        print(contents, file=sys.stderr)

    fileop = "appended" if append else "written"
    message = "{0} {1} to `{2}`.".format(meta, fileop, filename)
    logging.debug(message.capitalize())
    if meta == "run script" and not append:
        sh("chmod u+x {0}".format(filename))


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


def read_block(handle, signal):
    """
    Useful for reading block-like file formats, for example FASTA or OBO file,
    such file usually startswith some signal, and in-between the signals are a
    record
    """
    signal_len = len(signal)
    it = (
        x[1]
        for x in groupby(handle, key=lambda row: row.strip()[:signal_len] == signal)
    )
    found_signal = False
    for header in it:
        header = list(header)
        for h in header[:-1]:
            h = h.strip()
            if h[:signal_len] != signal:
                continue
            yield h, []  # Header only, no contents
        header = header[-1].strip()
        if header[:signal_len] != signal:
            continue
        found_signal = True
        seq = list(s.strip() for s in next(it))
        yield header, seq

    if not found_signal:
        handle.seek(0)
        seq = list(s.strip() for s in handle)
        yield None, seq


def is_number(s, cast=float):
    """
    Check if a string is a number. Use cast=int to check if s is an integer.
    """
    try:
        cast(s)  # for int, long and float
    except ValueError:
        return False

    return True


def get_number(s, cast=int):
    """
    Try to get a number out of a string, and cast it.
    """
    import string

    d = "".join(x for x in str(s) if x in string.digits)
    return cast(d) if d else s


def flexible_cast(s):
    if is_number(s, cast=int):
        return int(s)
    elif is_number(s, cast=float):
        return float(s)
    return s


def main():

    actions = (
        ("pairwise", "convert a list of IDs into all pairs"),
        ("split", "split large file into N chunks"),
        ("reorder", "reorder columns in tab-delimited files"),
        ("flatten", "convert a list of IDs into one per line"),
        ("unflatten", "convert lines to a list of IDs on single line"),
        ("group", "group elements in a table based on key (groupby) column"),
        ("setop", "set operations on files"),
        ("join", "join tabular-like files based on common column"),
        ("subset", "subset tabular-like files based on common column"),
        ("truncate", "remove lines from end of file"),
        ("append", "append a column with fixed value"),
        ("seqids", "make a list of seqids for graphics.karyotype"),
        ("mergecsv", "merge a set of tsv files"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def seqids(args):
    """
    %prog seqids prefix start end

    Make a list of seqids for graphics.karyotype. For example:

    $ python -m jcvi.formats.base seqids chromosome_ 1 3
    chromosome_1,chromosome_2,chromosome_3
    $ python -m jcvi.formats.base seqids A 3 1 --pad0=2
    A03,A02,A01
    """
    p = OptionParser(seqids.__doc__)
    p.add_option("--pad0", default=0, help="How many zeros to pad")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    prefix, start, end = args
    pad0 = opts.pad0
    start, end = int(start), int(end)
    step = 1 if start <= end else -1

    print(
        ",".join(
            [
                "{}{:0{}d}".format(prefix, x, pad0)
                for x in range(start, end + step, step)
            ]
        )
    )


def pairwise(args):
    """
    %prog pairwise ids

    Convert a list of IDs into all pairs.
    """
    from itertools import combinations

    p = OptionParser(pairwise.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (idsfile,) = args
    ids = SetFile(idsfile)
    ids = sorted(ids)
    fw = open(idsfile + ".pairs", "w")
    for a, b in combinations(ids, 2):
        print("\t".join((a, b)), file=fw)
    fw.close()


def append(args):
    """
    %prog append csvfile [tag]

    Append a column with fixed value. If tag is missing then just append the
    filename.
    """
    p = OptionParser(append.__doc__)
    p.set_sep()
    p.set_outfile()
    opts, args = p.parse_args(args)

    nargs = len(args)
    if nargs not in (1, 2):
        sys.exit(not p.print_help())

    csvfile = args[0]
    tag = args[1] if nargs == 2 else csvfile
    fp = must_open(csvfile)
    fw = must_open(opts.outfile, "w")
    for row in fp:
        row = row.rstrip("\r\n")
        row = opts.sep.join((row, tag))
        print(row, file=fw)


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
    while f.tell() > 0:
        f.seek(-1, os.SEEK_CUR)
        char = f.read(1)
        if char == "\n":
            count += 1
        if count == number + 1:
            f.truncate()
            print("Removed {0} lines from end of file".format(number), file=sys.stderr)
            return number

        f.seek(-1, os.SEEK_CUR)

    if count < number + 1:
        print("No change: requested removal would leave empty file", file=sys.stderr)
        return -1


def flatten(args):
    """
    %prog flatten filename > ids

    Convert a list of IDs (say, multiple IDs per line) and move them into one
    per line.

    For example, convert this, to this:
    A,B,C                    | A
    1                        | B
    a,4                      | C
                             | 1
                             | a
                             | 4

    If multi-column file with multiple elements per column, zip then flatten like so:
    A,B,C    2,10,gg         | A,2
    1,3      4               | B,10
                             | C,gg
                             | 1,4
                             | 3,na
    """
    from itertools import zip_longest

    p = OptionParser(flatten.__doc__)
    p.set_sep(sep=",")
    p.add_option(
        "--zipflatten",
        default=None,
        dest="zipsep",
        help="Specify if columns of the file should be zipped before"
        + " flattening. If so, specify delimiter separating column elements",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (tabfile,) = args
    zipsep = opts.zipsep

    fp = must_open(tabfile)
    for row in fp:
        if zipsep:
            row = row.rstrip()
            atoms = row.split(opts.sep)
            frows = []
            for atom in atoms:
                frows.append(atom.split(zipsep))
            print(
                "\n".join(
                    [zipsep.join(x) for x in list(zip_longest(*frows, fillvalue="na"))]
                )
            )
        else:
            print(row.strip().replace(opts.sep, "\n"))


def unflatten(args):
    """
    %prog unflatten idsfile > unflattened

    Given a list of ids, one per line, unflatten the list onto a single line with sep.
    """
    p = OptionParser(unflatten.__doc__)
    p.add_option("--sep", default=",", help="Separator when joining ids")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (idsfile,) = args
    ids = must_open(idsfile).read().split()
    with must_open(opts.outfile, "w") as fw:
        print(opts.sep.join(ids), file=fw)


def group(args):
    """
    %prog group tabfile > tabfile.grouped

    Given a tab-delimited file, either group all elements within the file or
    group the elements in the value column(s) based on the key (groupby) column

    For example, convert this | into this
    ---------------------------------------
    a   2    3    4           | a,2,3,4,5,6
    a   5    6                | b,7,8
    b   7    8                | c,9,10,11
    c   9                     |
    c  10   11                |

    If grouping by a particular column,
    convert this              | into this:
    ---------------------------------------------
    a   2    3    4           | a   2,5   3,6   4
    a   5    6                | b   7     8
    b   7    8                | c   9,10  11
    c   9                     |
    c  10   11                |

    By default, it uniqifies all the grouped elements
    """
    from jcvi.utils.cbook import AutoVivification
    from jcvi.utils.grouper import Grouper

    p = OptionParser(group.__doc__)
    p.set_sep()
    p.add_option(
        "--groupby", default=None, type="int", help="Default column to groupby"
    )
    p.add_option(
        "--groupsep", default=",", help="Separator to join the grouped elements"
    )
    p.add_option(
        "--nouniq",
        default=False,
        action="store_true",
        help="Do not uniqify the grouped elements",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (tabfile,) = args
    sep = opts.sep
    groupby = opts.groupby
    groupsep = opts.groupsep

    cols = []
    grouper = AutoVivification() if groupby is not None else Grouper()
    fp = must_open(tabfile)
    for row in fp:
        row = row.rstrip()
        atoms = row.split(sep)
        if groupby is not None:
            if len(cols) < len(atoms):
                cols = [x for x in range(len(atoms))]
            if groupby not in cols:
                logging.error("groupby col index `{0}` is out of range".format(groupby))
                sys.exit()

            key = atoms[groupby]
            for col in cols:
                if col == groupby:
                    continue
                if not grouper[key][col]:
                    grouper[key][col] = [] if opts.nouniq else set()
                if col < len(atoms):
                    if groupsep in atoms[col]:
                        for atom in atoms[col].split(groupsep):
                            if opts.nouniq:
                                grouper[key][col].append(atom)
                            else:
                                grouper[key][col].add(atom)
                    else:
                        if opts.nouniq:
                            grouper[key][col].append(atoms[col])
                        else:
                            grouper[key][col].add(atoms[col])
        else:
            grouper.join(*atoms)

    for key in grouper:
        if groupby is not None:
            line = []
            for col in cols:
                if col == groupby:
                    line.append(key)
                elif col in grouper[key].keys():
                    line.append(groupsep.join(grouper[key][col]))
                else:
                    line.append("na")
            print(sep.join(line))
        else:
            print(groupsep.join(key))


def reorder(args):
    """
    %prog reorder tabfile 1,2,4,3 > newtabfile

    Reorder columns in tab-delimited files. The above syntax will print out a
    new file with col-1,2,4,3 from the old file.
    """
    import csv

    p = OptionParser(reorder.__doc__)
    p.set_sep()
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
    %prog split file outdir N

    Split file into N records. This allows splitting FASTA/FASTQ/TXT file
    properly at boundary of records. Split is useful for parallelization
    on input chunks.

    Option --mode is useful on how to break into chunks.
    1. chunk - chunk records sequentially, 1-100 in file 1, 101-200 in file 2, etc.
    2. cycle - chunk records in Round Robin fashion
    3. optimal - try to make split file of roughly similar sizes, using LPT
    algorithm. This is the default.
    """
    p = OptionParser(split.__doc__)
    mode_choices = ("batch", "cycle", "optimal")
    p.add_option("--all", default=False, action="store_true", help="split all records")
    p.add_option(
        "--mode",
        default="optimal",
        choices=mode_choices,
        help="Mode when splitting records",
    )
    p.add_option(
        "--format", choices=("fasta", "fastq", "txt", "clust"), help="input file format"
    )

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    filename, outdir, N = args
    fs = FileSplitter(filename, outputdir=outdir, format=opts.format, mode=opts.mode)

    if opts.all:
        logging.debug("option -all override N")
        N = fs.num_records
    else:
        N = min(fs.num_records, int(N))
        assert N > 0, "N must be > 0"

    logging.debug("split file into %d chunks" % N)
    fs.split(N)

    return fs


def join(args):
    """
    %prog join file1.txt(pivotfile) file2.txt ..

    Join tabular-like files based on common column.
    --column specifies the column index to pivot on.
      Use comma to separate multiple values if the pivot column is different
      in each file. Maintain the order in the first file.
    --sep specifies the column separators, default to tab.
      Use comma to separate multiple values if the column separator is different
      in each file.
    """
    p = OptionParser(join.__doc__)
    p.add_option(
        "--column", default="0", help="0-based column id, multiple values allowed"
    )
    p.set_sep(multiple=True)
    p.add_option(
        "--noheader", default=False, action="store_true", help="Do not print header"
    )
    p.add_option("--na", default="na", help="Value for unjoined data")
    p.add_option(
        "--compact",
        default=False,
        action="store_true",
        help="Do not repeat pivotal columns in output",
    )
    p.add_option(
        "--keysep",
        default=",",
        help="specify separator joining multiple elements in the key column"
        + " of the pivot file",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)
    nargs = len(args)

    keysep = opts.keysep
    compact = opts.compact

    if len(args) < 2:
        sys.exit(not p.print_help())

    na = opts.na
    c = opts.column
    if "," in c:
        cc = [int(x) for x in c.split(",")]
    else:
        cc = [int(c)] * nargs

    assert len(cc) == nargs, "Column index number != File number"

    s = opts.sep
    if "," in s:
        ss = [x for x in s.split(",")]
    else:
        ss = [s] * nargs

    assert len(ss) == nargs, "column separator number != File number"

    # Maintain the first file line order, and combine other files into it
    pivotfile = args[0]
    files = [
        DictFile(f, keypos=c, valuepos=None, delimiter=s)
        for f, c, s in zip(args, cc, ss)
    ]
    otherfiles = files[1:]
    # The header contains filenames
    headers = []
    for i, x in enumerate(files):
        ncols = x.ncols
        if i and compact:
            ncols -= 1
        headers += [op.basename(x.filename)] * ncols
    header = "\t".join(headers)

    fp = must_open(pivotfile)
    fw = must_open(opts.outfile, "w")
    if not opts.noheader:
        print(header, file=fw)

    for row in fp:
        row = row.rstrip()
        atoms = row.split(ss[0])
        newrow = atoms
        key = atoms[cc[0]]
        keys = key.split(keysep) if keysep in key else [key]
        for d in otherfiles:
            drows = list()
            for key in keys:
                krow = d.get(key, [na] * d.ncols)
                if compact:
                    krow.pop(d.keypos)
                drows.append(krow)
            drow = [keysep.join(x) for x in list(zip(*drows))]
            newrow += drow
        print("\t".join(newrow), file=fw)


def subset(args):
    """
    %prog subset file1.txt(pivotfile) file2.txt ..

    subset tabular-like file1 based on common column with file 2.
    Normally file1 should have unique row entries.
    If more than one file2 are provided, they must have same column separators.
    Multiple file2's will be concatenated in the output.

    --column specifies the column index (0-based) to pivot on.
      Use comma to separate multiple values if the pivot column is different
      in each file. Maintain the order in the first file.
    --sep specifies the column separators, default to tab.
      Use comma to separate multiple values if the column separator is different
      in each file.
    """

    p = OptionParser(subset.__doc__)
    p.add_option(
        "--column", default="0", help="0-based column id, multiple values allowed"
    )
    p.set_sep(multiple=True)
    p.add_option(
        "--pivot",
        default=1,
        type="int",
        help="1 for using order in file1, 2 for using order in \
                    file2",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)
    nargs = len(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    c = opts.column
    if "," in c:
        cc = [int(x) for x in c.split(",")]
        assert len(set(cc[1:])) == 1, "Multiple file2's must have same column index."
        cc = cc[0:2]
    else:
        cc = [int(c)] * 2

    s = opts.sep
    if "," in s:
        ss = [x for x in s.split(",")]
        assert (
            len(set(cc[1:])) == 1
        ), "Multiple file2's must have same column separator."
        ss = ss[0:2]
    else:
        ss = [s] * 2

    if nargs > 2:
        file2 = FileMerger(args[1:], outfile="concatenatedFile2").merge()
    else:
        file2 = args[1]
    newargs = [args[0], file2]

    files = [
        DictFile(f, keypos=c, valuepos=None, delimiter=s)
        for f, c, s in zip(newargs, cc, ss)
    ]

    pivot = 0 if opts.pivot == 1 else 1
    fp = open(newargs[pivot])
    fw = must_open(opts.outfile, "w")

    for row in fp:
        row = row.rstrip()
        atoms = row.split(ss[pivot])
        key = atoms[cc[pivot]]
        d = files[1 - pivot]
        if key in d:
            print(ss[0].join(files[0][key]), file=fw)

    if nargs > 2:
        FileShredder([file2])


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
    from natsort import natsorted

    p = OptionParser(setop.__doc__)
    p.add_option(
        "--column",
        default=0,
        type="int",
        help="The column to extract, 0-based, -1 to disable",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (statement,) = args
    fa, op, fb = statement.split()
    assert op in ("|", "&", "-", "^")

    column = opts.column
    fa = SetFile(fa, column=column)
    fb = SetFile(fb, column=column)

    if op == "|":
        t = fa | fb
    elif op == "&":
        t = fa & fb
    elif op == "-":
        t = fa - fb
    elif op == "^":
        t = fa ^ fb

    for x in natsorted(t):
        print(x)


def mergecsv(args):
    """
    %prog mergecsv *.tsv

    Merge a set of tsv files.
    """
    p = OptionParser(mergecsv.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    tsvfiles = args
    outfile = opts.outfile

    if op.exists(outfile):
        os.remove(outfile)

    tsvfile = tsvfiles[0]
    fw = must_open(opts.outfile, "w")
    for i, tsvfile in enumerate(tsvfiles):
        fp = open(tsvfile)
        if i > 0:
            next(fp)
        for row in fp:
            fw.write(row)
    fw.close()


if __name__ == "__main__":
    main()
