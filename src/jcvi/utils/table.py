"""
Routines to summarize and report tabular data.
"""


def comment_banner(s, width=50):
    line = "#" * width
    return "\n".join((line, "#", "# " + s.strip(), "#", line))


def banner(header, rows, major="=", minor="-"):
    formatted = [header] + rows
    rulersize = max(max(len(z) for z in x.splitlines()) for x in formatted)
    table_edge = major * rulersize
    table_sep = minor * rulersize
    rows = "\n".join(rows)

    return "\n".join((table_edge, header, table_sep, rows, table_sep))


def loadtable(header, rows, thousands=True):
    """
    Print a tabular output, with horizontal separators
    """
    formatted = load_csv(header, rows, sep="   ", thousands=thousands)
    header, rows = formatted[0], formatted[1:]

    return banner(header, rows)


def tabulate(d, transpose=False, thousands=True, key_fun=None, sep=",", align=True):
    """
    d is a dictionary, keyed by tuple(A, B).
    Goal is to put A in rows, B in columns, report data in table form.

    >>> d = {(1,'a'):3, (1,'b'):4, (2,'a'):5, (2,'b'):0}
    >>> print(tabulate(d))
    ===========
    o    a    b
    -----------
    1    3    4
    2    5    0
    -----------
    >>> print(tabulate(d, transpose=True))
    ===========
    o    1    2
    -----------
    a    3    5
    b    4    0
    -----------
    """
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
        data = [d.get(x, "n/a") for x in combo]
        data = ["{0:.1f}".format(x) if isinstance(x, float) else x for x in data]
        if key_fun:
            data = [key_fun(x) for x in data]
        table.append([str(r)] + data)

    if not align:
        formatted = load_csv(header, table, sep=sep)
        return "\n".join(formatted)

    return loadtable(header, table, thousands=thousands)


def load_csv(header, contents, sep=",", thousands=False, align=True):

    from jcvi.formats.base import is_number
    from jcvi.utils.cbook import thousands as th

    allcontents = [header] + contents if header else contents
    cols = len(contents[0])
    for content in allcontents:
        assert len(content) == cols

    # Stringify the contents
    for i, content in enumerate(allcontents):
        if thousands:
            content = [int(x) if is_number(x, cast=int) else x for x in content]
            content = [
                th(x) if (is_number(x, cast=int) and x >= 1000) else x for x in content
            ]
        allcontents[i] = [str(x) for x in content]

    colwidths = [max(len(x[i]) for x in allcontents) for i in range(cols)]
    sep += " "
    formatted_contents = []
    for content in allcontents:
        rjusted = (
            [x.rjust(cw) for x, cw in zip(content, colwidths)] if align else content
        )
        formatted = sep.join(rjusted)
        formatted_contents.append(formatted)

    return formatted_contents


def write_csv(
    header,
    contents,
    sep=",",
    filename="stdout",
    thousands=False,
    tee=False,
    align=True,
    comment=False,
):
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

    formatted = load_csv(header, contents, sep=sep, thousands=thousands, align=align)
    if comment:
        formatted[0] = "#" + formatted[0][1:]
    formatted = "\n".join(formatted)
    output = must_open(filename, "w")
    print(formatted, file=output)
    if tee and filename != "stdout":
        print(formatted)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
