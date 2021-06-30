#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Manipulate PDF files, using PyPDF2 library.
"""
import os
import sys
import logging
import traceback

from natsort import natsorted

from PyPDF2 import PdfFileMerger, parse_filename_page_ranges
from PyPDF2.pagerange import PAGE_RANGE_HELP
from jcvi.formats.base import must_open
from jcvi.apps.base import OptionParser, ActionDispatcher


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
    p.add_option(
        "--nosort", default=False, action="store_true", help="Do not sort file names"
    )
    p.add_option(
        "--cleanup",
        default=False,
        action="store_true",
        help="Remove individual pdfs after merging",
    )
    p.set_outfile()
    p.set_verbose(help="Show page ranges as they are being read")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    outfile = opts.outfile
    if outfile in args:
        args.remove(outfile)

    if not opts.nosort:
        args = natsorted(args)

    filename_page_ranges = parse_filename_page_ranges(args)
    verbose = opts.verbose
    fw = must_open(outfile, "wb")

    merger = PdfFileMerger()
    in_fs = {}
    try:
        for (filename, page_range) in filename_page_ranges:
            if verbose:
                print(filename, page_range, file=sys.stderr)
            if filename not in in_fs:
                in_fs[filename] = open(filename, "rb")
            merger.append(in_fs[filename], pages=page_range)
    except:
        print(traceback.format_exc(), file=sys.stderr)
        print("Error while reading " + filename, file=sys.stderr)
        sys.exit(1)
    merger.write(fw)
    fw.close()

    if opts.cleanup:
        logging.debug("Cleaning up {} files".format(len(args)))
        for arg in args:
            os.remove(arg)


if __name__ == "__main__":
    main()
