#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Parse html pages.
"""
import os.path as op
import sys
import logging

from BeautifulSoup import BeautifulSoup
from urllib.parse import urljoin

from jcvi.apps.base import OptionParser, ActionDispatcher, download


def main():

    actions = (
        ("table", "convert HTML tables to csv"),
        ("links", "extract all links from web page"),
        ("gallery", "convert a folder of figures to a HTML table"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def gallery(args):
    """
    %prog gallery folder link_prefix

    Convert a folder of figures to a HTML table. For example:

    $ python -m jcvi.formats.html gallery Paper-figures/
    https://dl.dropboxusercontent.com/u/15937715/Data/Paper-figures/

    Maps the images from local to remote.
    """
    from more_itertools import grouper
    from jcvi.apps.base import iglob

    p = OptionParser(gallery.__doc__)
    p.add_option("--columns", default=3, type="int", help="How many cells per row")
    p.add_option("--width", default=200, type="int", help="Image width")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, link_prefix = args
    width = opts.width
    images = iglob(folder, "*.jpg,*.JPG,*.png")
    td = '<td>{0}<br><a href="{1}"><img src="{1}" width="{2}"></a></td>'
    print("<table>")
    for ims in grouper(images, opts.columns):
        print('<tr height="{0}" valign="top">'.format(width + 5))
        for im in ims:
            if not im:
                continue
            im = op.basename(im)
            pf = im.split(".")[0].replace("_", "-")
            link = link_prefix.rstrip("/") + "/" + im
            print(td.format(pf, link, width))
        print("</tr>")
    print("</table>")


def links(args):
    """
    %prog links url

    Extract all the links "<a href=''>" from web page.
    """
    p = OptionParser(links.__doc__)
    p.add_option(
        "--img",
        default=False,
        action="store_true",
        help="Extract <img> tags [default: %default]",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (url,) = args
    img = opts.img

    htmlfile = download(url)
    page = open(htmlfile).read()
    soup = BeautifulSoup(page)

    tag = "img" if img else "a"
    src = "src" if img else "href"
    aa = soup.findAll(tag)
    for a in aa:
        link = a.get(src)
        link = urljoin(url, link)
        print(link)


def unescape(s, unicode_action="replace"):
    """
    Unescape HTML strings, and convert &amp; etc.
    """
    from html.parser import HTMLParser

    hp = HTMLParser.HTMLParser()
    s = hp.unescape(s)
    s = s.encode("ascii", unicode_action)
    s = s.replace("\n", "").strip()
    return s


def table(args):
    """
    %prog table page.html

    Convert HTML tables to csv.
    """
    import csv

    p = OptionParser(table.__doc__)
    p.set_sep(sep=",")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (htmlfile,) = args
    page = open(htmlfile).read()
    soup = BeautifulSoup(page)

    for i, tabl in enumerate(soup.findAll("table")):
        nrows = 0
        csvfile = htmlfile.rsplit(".", 1)[0] + ".{0}.csv".format(i)
        writer = csv.writer(open(csvfile, "w"), delimiter=opts.sep)
        rows = tabl.findAll("tr")
        for tr in rows:
            cols = tr.findAll("td")
            if not cols:
                cols = tr.findAll("th")

            row = []
            for td in cols:
                try:
                    cell = "".join(td.find(text=True))
                    cell = unescape(cell)
                except TypeError:
                    cell = ""
                row.append(cell)
            writer.writerow(row)
            nrows += 1
        logging.debug("Table with {0} rows written to `{1}`.".format(nrows, csvfile))


if __name__ == "__main__":
    main()
