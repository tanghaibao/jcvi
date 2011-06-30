"""
Routines to summarize and report tabular data.
"""


def banner(listOfStuff, rulersize=80, major='=', minor='-'):
    """
    Print a tabular output, with horizontal separators
    """
    table_edge = major * rulersize + "\n"
    table_sep = minor * rulersize + "\n"
    contents = table_sep.join(str(x) + "\n" for x in listOfStuff)
    return "".join((table_edge, contents, table_edge))


def loadtable(header, rows):

    from cogent import LoadTable
    return LoadTable(header=header, row=row)


def tabulate(d, transpose=False, key_fun=None):
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
    >>> print tabulate(d, transpose=True)
    ===========
    o    1    2
    -----------
    a    3    5
    b    4    0
    -----------
    """
    from cogent import LoadTable

    pairs = d.keys()
    rows, cols = zip(*pairs)
    if transpose:
        rows, cols = cols, rows

    rows = sorted(set(rows))
    cols = sorted(set(cols))
    header = ["o"] + list(cols)
    table = []
    for r in rows:
        combo = [(r, c) for c in cols]
        if transpose:
            combo = [(c, r) for (r, c) in combo]
        data = [d[x] for x in combo]
        if key_fun:
            data = [key_fun(x) for x in data]
        table.append([str(r)] + data)

    table = LoadTable(header=header, rows=table)

    return table


if __name__ == '__main__':
    import doctest
    doctest.testmod()
