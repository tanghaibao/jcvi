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
    return LoadTable(header=header, rows=rows)


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
        data = ["{0:.1f}".format(x) if isinstance(x, float) else x \
                    for x in data]
        if key_fun:
            data = [key_fun(x) for x in data]
        table.append([str(r)] + data)

    table = LoadTable(header=header, rows=table)

    return table


def write_csv(header, contents, sep=",", filename="stdout", tee=False):
    """
    Write csv that are aligned with the column headers.

    >>> header = ["x_value", "y_value"]
    >>> contents = [(1, 100), (2, 200)]
    >>> write_csv(header, contents)
    x_value, y_value
          1,     100
          2,     200
    """
    from jcvi.formats.base import must_open

    fw = must_open(filename, "w")
    allcontents = [header] + contents if header else contents
    cols = len(contents[0])
    for content in allcontents:
        assert len(content) == cols

    # Stringify the contents
    for i, content in enumerate(allcontents):
        allcontents[i] = [str(x) for x in content]

    colwidths = [max(len(x[i]) for x in allcontents) for i in xrange(cols)]
    sep += " "
    for content in allcontents:
        rjusted = [x.rjust(cw) for x, cw in zip(content, colwidths)]
        formatted = sep.join(rjusted)
        print >> fw, formatted
        if tee and filename != "stdout":
            print formatted


if __name__ == '__main__':
    import doctest
    doctest.testmod()
