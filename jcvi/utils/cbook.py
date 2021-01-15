"""
Useful recipes from various internet sources (thanks)
mostly decorator patterns
"""
from __future__ import print_function

import os.path as op
import re
import sys
import logging
import functools

from collections import defaultdict


def inspect(item, maxchar=80):
    """
    Inspect the attributes of an item.
    """
    for i in dir(item):
        try:
            member = str(getattr(item, i))
            if maxchar and len(member) > maxchar:
                member = member[:maxchar] + "..."
        except:
            member = "[ERROR]"
        print("{}: {}".format(i, member), file=sys.stderr)


def timeit(func):
    """
    <http://www.zopyx.com/blog/a-python-decorator-for-measuring-the-execution-time-of-methods>
    """
    import time

    def timed(*args, **kw):
        ts = time.time()
        result = func(*args, **kw)
        te = time.time()

        msg = "{0}{1} {2:.2f}s".format(func.__name__, args, te - ts)
        logging.debug(msg)

        return result

    return timed


def depends(func):
    """
    Decorator to perform check on infile and outfile. When infile is not present, issue
    warning, and when outfile is present, skip function calls.
    """
    from jcvi.apps.base import need_update, listify

    infile = "infile"
    outfile = "outfile"

    def wrapper(*args, **kwargs):
        assert outfile in kwargs, "You need to specify `outfile=` on function call"
        if infile in kwargs:
            infilename = listify(kwargs[infile])
            for x in infilename:
                assert op.exists(x), "The specified infile `{0}` does not exist".format(
                    x
                )

        outfilename = kwargs[outfile]
        if need_update(infilename, outfilename):
            return func(*args, **kwargs)
        else:
            msg = "File `{0}` exists. Computation skipped.".format(outfilename)
            logging.debug(msg)

        outfilename = listify(outfilename)

        for x in outfilename:
            assert op.exists(x), "Something went wrong, `{0}` not found".format(x)

        return outfilename

    return wrapper


"""
Functions that make text formatting easier.
"""


class Registry(defaultdict):
    def __init__(self, *args, **kwargs):
        super(Registry, self).__init__(list, *args, **kwargs)

    def iter_tag(self, tag):
        for key, ts in self.items():
            if tag in ts:
                yield key

    def get_tag(self, tag):
        return list(self.iter_tag(tag))

    def count(self, tag):
        return sum(1 for x in self.iter_tag(tag))

    def update_from(self, filename):
        from jcvi.formats.base import DictFile

        d = DictFile(filename)
        for k, v in d.items():
            self[k].append(v)


class SummaryStats(object):
    def __init__(self, a, dtype=None, title=None):
        import numpy as np

        self.data = a = np.array(a, dtype=dtype)
        self.min = a.min()
        self.max = a.max()
        self.size = a.size
        self.mean = np.mean(a)
        self.sd = np.std(a)
        self.median = np.median(a)
        self.sum = a.sum()
        self.title = title

        a.sort()
        self.firstq = a[self.size // 4]
        self.thirdq = a[self.size * 3 // 4]
        self.p1 = a[int(self.size * 0.025)]
        self.p2 = a[int(self.size * 0.975)]

        if dtype == "int":
            self.mean = int(self.mean)
            self.sd = int(self.sd)
            self.median = int(self.median)

    def __str__(self):
        s = self.title + ": " if self.title else ""
        s += "Min={} Max={} N={} Mean={:.2f} SD={:.2f} Median={} Sum={}".format(
            self.min, self.max, self.size, self.mean, self.sd, self.median, self.sum
        )
        return s

    def todict(self, quartile=False):
        d = {"Min": self.min, "Max": self.max, "Mean": self.mean, "Median": self.median}
        if quartile:
            d.update({"1st Quartile": self.firstq, "3rd Quartile": self.thirdq})

        return d

    def tofile(self, filename):
        fw = open(filename, "w")
        for x in self.data:
            print(x, file=fw)
        fw.close()
        logging.debug(
            "Array of size {0} written to file `{1}`.".format(self.size, filename)
        )


class AutoVivification(dict):
    """
    Implementation of perl's autovivification feature.

    Thanks to <http://stackoverflow.com/questions/651794/whats-the-best-way-to-initialize-a-dict-of-dicts-in-python>
    """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def enumerate_reversed(sequence):
    """
    Perform reverse enumeration, returning an iterator with decrementing
    index/position values

    Source: http://stackoverflow.com/questions/529424/traverse-a-list-in-reverse-order-in-python
    """
    for index in reversed(range(len(sequence))):
        yield index, sequence[index]


def percentage(a, b, precision=1, mode=0):
    """
    >>> percentage(100, 200)
    '100 of 200 (50.0%)'
    """
    _a, _b = a, b
    pct = "{0:.{1}f}%".format(a * 100.0 / b, precision)
    a, b = thousands(a), thousands(b)
    if mode == 0:
        return "{0} of {1} ({2})".format(a, b, pct)
    elif mode == 1:
        return "{0} ({1})".format(a, pct)
    elif mode == 2:
        return _a * 100.0 / _b
    return pct


def thousands(x):
    """
    >>> thousands(12345)
    '12,345'
    """
    import locale

    try:
        locale.setlocale(locale.LC_ALL, "en_US.utf8")
    except Exception:
        locale.setlocale(locale.LC_ALL, "en_US.UTF-8")
    finally:
        s = "%d" % x
        groups = []
        while s and s[-1].isdigit():
            groups.append(s[-3:])
            s = s[:-3]
        return s + ",".join(reversed(groups))
    return locale.format("%d", x, True)


SUFFIXES = {
    1000: ["", "Kb", "Mb", "Gb", "Tb", "Pb", "Eb", "Zb"],
    1024: ["B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB"],
}


def human_size(size, a_kilobyte_is_1024_bytes=False, precision=1, target=None):
    """Convert a file size to human-readable form.

    Keyword arguments:
    size -- file size in bytes
    a_kilobyte_is_1024_bytes -- if True (default), use multiples of 1024
                                if False, use multiples of 1000

    Returns: string
    Credit: <http://diveintopython3.org/your-first-python-program.html>

    >>> print(human_size(1000000000000, True))
    931.3GiB
    >>> print(human_size(1000000000000))
    1.0Tb
    >>> print(human_size(300))
    300.0
    """
    if size < 0:
        raise ValueError("number must be non-negative")

    multiple = 1024 if a_kilobyte_is_1024_bytes else 1000
    for suffix in SUFFIXES[multiple]:

        if target:
            if suffix == target:
                break
            size /= float(multiple)
        else:
            if size >= multiple:
                size /= float(multiple)
            else:
                break

    return "{0:.{1}f}{2}".format(size, precision, suffix)


def autoscale(bp, optimal=6):
    """
    >>> autoscale(150000000)
    20000000
    >>> autoscale(97352632)
    10000000
    """
    slen = str(bp)
    tlen = slen[0:2] if len(slen) > 1 else slen[0]
    precision = len(slen) - 2  # how many zeros we need to pad?
    bp_len_scaled = int(tlen)  # scale bp_len to range (0, 100)
    tick_diffs = [(x, abs(bp_len_scaled / x - optimal)) for x in [1, 2, 5, 10]]
    best_stride, best_tick_diff = min(tick_diffs, key=lambda x: x[1])

    while precision > 0:
        best_stride *= 10
        precision -= 1

    return best_stride


def gene_name(st, exclude=("ev",), sep="."):
    """
    Helper functions in the BLAST filtering to get rid alternative splicings.
    This is ugly, but different annotation groups are inconsistent with respect
    to how the alternative splicings are named. Mostly it can be done by removing
    the suffix, except for ones in the exclude list.
    """
    if any(st.startswith(x) for x in exclude):
        sep = None
    st = st.split("|")[0]

    if sep and sep in st:
        name, suffix = st.rsplit(sep, 1)
    else:
        name, suffix = st, ""

    # We only want to remove suffix that are isoforms, longer suffix would
    # suggest that it is part of the right gene name
    if len(suffix) != 1:
        name = st

    return name


def seqid_parse(seqid, sep=["-"], stdpf=True):
    """
    This function tries to parse seqid (1st col in bed files)
    return prefix, numeric id, and suffix, for example:

    >>> seqid_parse('chr1_random')
    ('Chr', '1', '_random')
    >>> seqid_parse('AmTr_v1.0_scaffold00001', '', stdpf=False)
    ('AmTr_v1.0_scaffold', '00001', '')
    >>> seqid_parse('AmTr_v1.0_scaffold00001')
    ('Sca', '00001', '')
    >>> seqid_parse('PDK_30s1055861')
    ('C', '1055861', '')
    >>> seqid_parse('PDK_30s1055861', stdpf=False)
    ('PDK', '1055861', '')
    >>> seqid_parse("AC235758.1", stdpf=False)
    ('AC', '235758.1', '')
    """
    seqid = seqid.split(";")[0]
    if "mito" in seqid or "chloro" in seqid:
        return (seqid, "", "")

    numbers = re.findall(r"\d+\.*\d*", seqid)

    if not numbers:
        return (seqid, "", "")

    id = numbers[-1]
    lastnumi = seqid.rfind(id)
    suffixi = lastnumi + len(id)
    suffix = seqid[suffixi:]

    if sep is None:
        sep = [""]
    elif type(sep) == str:
        sep = [sep]

    prefix = seqid[:lastnumi]
    if not stdpf:
        sep = "|".join(sep)
        atoms = re.split(sep, prefix)
        if len(atoms) == 1:
            prefix = atoms[0]
        else:
            prefix = atoms[-2]
        prefix = prefix.replace("Chromosome", "Chr")
    else:  # use standard prefix
        if re.findall("chr", prefix, re.I):
            prefix = "Chr"
        if re.findall("lg", prefix, re.I):
            prefix = "LG"
        elif re.findall("sca", prefix, re.I):
            prefix = "Sca"
        elif re.findall("supercontig", prefix, re.I):
            prefix = "SCg"
        elif re.findall("ctg|contig", prefix, re.I):
            prefix = "Ctg"
        elif re.findall("BAC", prefix, re.I):
            prefix = "BAC"
        else:
            prefix = "C"

    return prefix, id, suffix


def fixChromName(name, orgn="medicago"):
    """
    Convert quirky chromosome names encountered in different
    release files, which are very project specific, into a more
    general format.

    For example, in Medicago
        Convert a seqid like
            `Mt3.5.1_Chr1` to `chr1`
            `Mt3.5_Chr3` to `chr3`
            `chr01_pseudomolecule_IMGAG` to `chr1`

    Some examples from Maize
        Convert a seqid like
            `chromosome:AGPv2:2:1:237068873:1` to `2`
        Special cases
            `chromosome:AGPv2:mitochondrion:1:569630:1` to `Mt`
            `chromosome:AGPv2:chloroplast:1:140384:1` to `Pt`
    """
    import re

    mtr_pat1 = re.compile(r"Mt[0-9]+\.[0-9]+[\.[0-9]+]{0,}_([a-z]+[0-9]+)")
    mtr_pat2 = re.compile(r"([A-z0-9]+)_[A-z]+_[A-z]+")

    zmays_pat = re.compile(r"[a-z]+:[A-z0-9]+:([A-z0-9]+):[0-9]+:[0-9]+:[0-9]+")
    zmays_sub = {"mitochondrion": "Mt", "chloroplast": "Pt"}
    if orgn == "medicago":
        for mtr_pat in (mtr_pat1, mtr_pat2):
            match = re.search(mtr_pat, name)
            if match:
                n = match.group(1)
                n = n.replace("0", "")
                name = re.sub(mtr_pat, n, name)
    elif orgn == "maize":
        match = re.search(zmays_pat, name)
        if match:
            n = match.group(1)
            name = re.sub(zmays_pat, n, name)
            if name in zmays_sub:
                name = zmays_sub[name]

    return name


def fill(text, delimiter="", width=70):
    """
    Wrap text with width per line
    """
    texts = []
    for i in range(0, len(text), width):
        t = delimiter.join(text[i : i + width])
        texts.append(t)
    return "\n".join(texts)


def tile(lt, width=70, gap=1):
    """
    Pretty print list of items.
    """
    from more_itertools import grouper

    max_len = max(len(x) for x in lt) + gap
    items_per_line = max(width // max_len, 1)
    lt = [x.rjust(max_len) for x in lt]
    g = list(grouper(lt, items_per_line, fillvalue=""))

    return "\n".join("".join(x) for x in g)


def uniqify(L):
    """
    Uniqify a list, maintains order (the first occurrence will be kept).
    """
    seen = set()
    nL = []
    for a in L:
        if a in seen:
            continue
        nL.append(a)
        seen.add(a)

    return nL


if __name__ == "__main__":
    import doctest

    doctest.testmod()
