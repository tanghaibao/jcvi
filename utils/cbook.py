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


def depends(func):
    """
    Decorator to perform check on infile and outfile. When infile is not present, issue
    warning, and when outfile is present, skip function calls.
    """
    from jcvi.apps.base import is_newer_file

    infile = "infile"
    outfile = "outfile"
    def wrapper(*args, **kwargs):
        assert outfile in kwargs, \
            "You need to specify `outfile=` on function call"
        if infile in kwargs:
            infilename = kwargs[infile]
            assert op.exists(infilename), \
                "The specified infile `{0}` does not exist" \
                    .format(infilename)

        outfilename = kwargs[outfile]
        if not op.exists(outfilename) or \
                is_newer_file(infilename, outfilename):
            func(*args, **kwargs)
        else:
            msg = "File `{0}` exists. Computation skipped." \
                .format(outfilename)
            logging.error(msg)

        assert op.exists(outfilename), \
                "Something went wrong, `{0}` not found" \
                .format(outfilename)

    return wrapper


"""
Functions that make text formatting easier.
"""

def percentage(a, b):
    """
    >>> percentage(100, 200)
    '100 of 200 (50.0%)'
    """
    return "{0} of {1} ({2:.1f}%)".format(a, b, a * 100. / b)


def thousands(x):
    """
    >>> thousands(12345)
    '12,345'
    """
    import locale
    locale.setlocale(locale.LC_ALL, "")
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
        if size >= multiple or (target and suffix != target):
            size /= float(multiple)
        else:
            return '{0:.{1}f}{2}'.format(size, precision, suffix)

    raise ValueError('number too large')


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
