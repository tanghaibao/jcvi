#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Manipulate PDF files, using PyPDF2 library.
"""
import sys

from natsort import natsorted

from PyPDF2 import PdfFileMerger, parse_filename_page_ranges
from PyPDF2.pagerange import PAGE_RANGE_HELP

from ..apps.base import ActionDispatcher, OptionParser, cleanup, logger

from .base import must_open


def main():

    actions = (("cat", "concatenate pages from pdf files into a single pdf file"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def cat(args):
    """
    %prog cat *.pdf -o output.pdf

    Concatenate pages from pdf files into a single pdf file.

    Page ranges refer to the previously-named file.
    A file not followed by a page range means all the pages of the file.

    PAGE RANGES are like Python slices.
            {page_range_help}
    EXAMPLES
        pdfcat -o output.pdf head.pdf content.pdf :6 7: tail.pdf -1
            Concatenate all of head.pdf, all but page seven of content.pdf,
            and the last page of tail.pdf, producing output.pdf.

        pdfcat chapter*.pdf >book.pdf
            You can specify the output file by redirection.

        pdfcat chapter?.pdf chapter10.pdf >book.pdf
            In case you don't want chapter 10 before chapter 2.
    """
    p = OptionParser(cat.__doc__.format(page_range_help=PAGE_RANGE_HELP))
    p.add_argument(
        "--nosort", default=False, action="store_true", help="Do not sort file names"
    )
    p.add_argument(
        "--cleanup",
        default=False,
        action="store_true",
        help="Remove individual pdfs after merging",
    )
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    outfile = opts.outfile
    if outfile in args:
        args.remove(outfile)

    should_sort = not opts.nosort
    if not all(x.endswith(".pdf") for x in args):
        should_sort = False
        logger.debug("Not sorting filenames because non-pdf args")

    if should_sort:
        args = natsorted(args)

    filename_page_ranges = parse_filename_page_ranges(args)
    nfiles = len(filename_page_ranges)
    merger = PdfFileMerger()
    with must_open(outfile, "wb") as fw:
        in_fs = {}
        try:
            for filename, page_range in filename_page_ranges:
                logger.debug("%s: %s", filename, page_range)
                if filename not in in_fs:
                    in_fs[filename] = open(filename, "rb")
                merger.append(in_fs[filename], pages=page_range)
        except Exception as e:
            logger.error("Error while reading %s: %s", filename, e)
            sys.exit(1)
        merger.write(fw)
        logger.info("Extracted %d files into `%s`", nfiles, outfile)

    if opts.cleanup:
        logger.debug("Cleaning up %d files", nfiles)
        cleanup(args)


if __name__ == "__main__":
    main()
