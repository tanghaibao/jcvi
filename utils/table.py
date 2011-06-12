"""
Routines to summarize and report tabular data.
"""

from cogent import LoadTable


def get_table():
    # http://foutaise.org/code/texttable/
    from texttable import Texttable
    table = Texttable()

    return table


def tabulate(d, key_fun=None):
    """
    d is a dictionary, keyed by tuple(A, B).
    Goal is to put A in rows, B in columns, report data in table form.

    >>> d = {(1,'a'):3, (1,'b'):4, (2,'a'):5, (2,'b'):0}
    >>> print tabulate(d)
    ===========
    o    a    b
    -----------
    1    3    4
    2    5    0
    -----------
    """
    pairs = d.keys()
    rows, cols = zip(*pairs)
    rows = sorted(set(rows))
    cols = sorted(set(cols))
    header = ["o"] + list(cols)
    table = []
    for r in rows:
        combo = [(r, c) for c in cols]
        data = [d[x] for x in combo]
        if key_fun:
            data = [key_fun(x) for x in data]
        table.append([str(r)] + data)

    table = LoadTable(header=header, rows=table)

    return table


if __name__ == '__main__':
    import doctest
    doctest.testmod()
