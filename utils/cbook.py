"""
Useful recipes from various internet sources (thanks)
mostly decorator patterns
"""

import os
import os.path as op
import logging


class memoized(object):
    """
    Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.

    Taken from recipe (http://wiki.python.org/moin/PythonDecoratorLibrary)
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        try:
            return self.cache[args]
        except KeyError:
            value = self.func(*args)
            self.cache[args] = value
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)

    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__

    def __get__(self, obj, objtype):
        """Support instance methods."""
        return functools.partial(self.__call__, obj)


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
    from jcvi.apps.base import need_update

    infile = "infile"
    outfile = "outfile"
    def wrapper(*args, **kwargs):
        assert outfile in kwargs, \
            "You need to specify `outfile=` on function call"
        if infile in kwargs:
            infilename = kwargs[infile]
            if isinstance(infilename, basestring):
                infilename = [infilename]
            for x in infilename:
                assert op.exists(x), \
                    "The specified infile `{0}` does not exist".format(x)

        outfilename = kwargs[outfile]
        if need_update(infilename, outfilename):
            return func(*args, **kwargs)
        else:
            msg = "File `{0}` exists. Computation skipped." \
                .format(outfilename)
            logging.debug(msg)

        if isinstance(outfilename, basestring):
            outfilename = [outfilename]

        for x in outfilename:
            assert op.exists(x), \
                    "Something went wrong, `{0}` not found".format(x)

        return outfilename

    return wrapper


"""
Functions that make text formatting easier.
"""

class SummaryStats (object):

    def __init__(self, a, title=None):
        import numpy as np

        self.data = a = np.array(a)
        self.min = a.min()
        self.max = a.max()
        self.size = a.size
        self.mean = np.mean(a)
        self.sd = np.std(a)
        self.median = np.median(a)
        self.title = title

        a.sort()
        self.firstq = a[self.size / 4]
        self.thirdq = a[self.size * 3 / 4]

    def __str__(self):
        s = self.title + ": " if self.title else ""
        s += "Min={0} Max={1} N={2} Mean={3:.0f} SD={4:.0f} Median={5:.0f}".\
                format(self.min, self.max, self.size,
                       self.mean, self.sd, self.median)
        return s

    def todict(self, quartile=False):
        d = {
            "Min": self.min, "Max": self.max,
            "Mean": self.mean, "Median": self.median
            }
        if quartile:
            d.update({
            "1st Quartile": self.firstq, "3rd Quartile": self.thirdq
            })

        return d

    def tofile(self, filename):
        fw = open(filename, "w")
        for x in self.data:
            print >> fw, x
        fw.close()
        logging.debug("Array of size {0} written to file `{1}`.".\
                        format(self.size, filename))


def percentage(a, b, denominator=True):
    """
    >>> percentage(100, 200)
    '100 of 200 (50.0%)'
    """
    if denominator:
        s = "{0} of {1} ({2:.1f}%)".format(a, b, a * 100. / b)
    else:
        s = "{0} ({1:.1f}%)".format(a, a * 100. / b)
    return s


def thousands(x):
    """
    >>> thousands(12345)
    '12,345'
    """
    import locale
    locale.setlocale(locale.LC_ALL, "en_AU.utf8")
    return locale.format('%d', x, True)


SUFFIXES = {1000: ['', 'Kb', 'Mb', 'Gb', 'Tb', 'Pb', 'Eb', 'Zb'],
            1024: ['B', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB', 'EiB', 'ZiB']}


def human_size(size, a_kilobyte_is_1024_bytes=False, precision=1, target=None):
    '''Convert a file size to human-readable form.

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
    '''
    if size < 0:
        raise ValueError('number must be non-negative')

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

    return '{0:.{1}f}{2}'.format(size, precision, suffix)


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

"""
Random ad-hoc functions
"""


def gene_name(st):
    """
    Helper functions in the BLAST filtering to get rid alternative splicings
    this is ugly, but different annotation groups are inconsistent
    with how the alternative splicings are named;
    mostly it can be done by removing the suffix
    except for papaya (evm...) and maize (somewhat complicated)
    """
    if st.startswith("ev"):
        return st
    if st.startswith("Os"):
        return st.rsplit("-", 1)[0]
    return st.rsplit(".", 1)[0]


def fill(text, delimiter="", width=70):
    """
    Wrap text with width per line
    """
    texts = []
    for i in xrange(0, len(text), width):
        t = delimiter.join(text[i:i + width])
        texts.append(t)
    return "\n".join(texts)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
