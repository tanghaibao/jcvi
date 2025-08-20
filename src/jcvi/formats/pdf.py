#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Manipulate PDF files, using PyPDF library.
"""
import shutil
import sys

from natsort import natsorted

from ..apps.base import ActionDispatcher, OptionParser, cleanup, logger
from .base import must_open

from pypdf import PdfWriter
from pypdf.pagerange import PageRange, parse_filename_page_ranges  # stable home


PAGE_RANGE_HELP = PageRange.__init__.__doc__


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

    # Fast-path: if exactly one input and range == ":" (all), do a simple copy
    if nfiles == 1:
        in_file, pr = filename_page_ranges[0]
        sl = getattr(pr, "to_slice", lambda: slice(None))()
        if sl.start is None and sl.stop is None and sl.step is None:
            logger.debug("Single input with full range detected; copying bytes")
            try:
                shutil.copyfile(in_file, outfile)
                logger.info("Copied `%s` to `%s`", in_file, outfile)
            except Exception as e:
                logger.error("Error while copying %s: %s", in_file, e)
                sys.exit(1)
            if opts.cleanup:
                logger.debug("Cleaning up 1 file")
                cleanup([in_file])
            return

    writer = PdfWriter()
    try:
        for filename, page_range in filename_page_ranges:
            logger.debug("%s: %s", filename, page_range)
            # `append` accepts PageRange objects directly (preferred over slices)
            writer.append(fileobj=filename, pages=page_range)

        with must_open(outfile, "wb") as fw:
            writer.write(fw)
        writer.close()
        logger.info("Extracted %d files into `%s`", nfiles, outfile)
    except Exception as e:
        logger.error("Error while reading %s: %s", filename, e)  # noqa: F821
        sys.exit(1)

    if opts.cleanup:
        logger.debug("Cleaning up %d files", nfiles)
        cleanup(args)


if __name__ == "__main__":
    main()
