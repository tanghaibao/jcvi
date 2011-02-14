"""
Useful recipes from various internet sources (thanks) 
mostly decorator patterns
"""

import logging


class memoized(object):
   """Decorator that caches a function's return value each time it is called.
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


SUFFIXES = {1000: ['Kb', 'Mb', 'Gb', 'Tb', 'Pb', 'Eb', 'Zb', 'Yb'],
            1024: ['KiB', 'MiB', 'GiB', 'TiB', 'PiB', 'EiB', 'ZiB', 'YiB']}


def human_size(size, a_kilobyte_is_1024_bytes=True, precision=1):
    '''Convert a file size to human-readable form.

    Keyword arguments:
    size -- file size in bytes
    a_kilobyte_is_1024_bytes -- if True (default), use multiples of 1024
                                if False, use multiples of 1000

    Returns: string
    Credit: <http://diveintopython3.org/your-first-python-program.html>
    '''
    if size < 0:
        raise ValueError('number must be non-negative')

    multiple = 1024 if a_kilobyte_is_1024_bytes else 1000
    for suffix in SUFFIXES[multiple]:
        size /= multiple
        if size < multiple:
            return '{0:.{1}f} {2}'.format(size, precision, suffix)

    raise ValueError('number too large')


if __name__ == '__main__':
    print(human_size(1000000000000, False))
    print(human_size(1000000000000))
