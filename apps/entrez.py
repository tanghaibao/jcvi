"""
Wrapper for calling Bio.Entrez tools to get the sequence from a list of IDs
"""

import os
import os.path as op
import sys
import time
import logging
import urllib2

from optparse import OptionParser
from Bio import Entrez, SeqIO

from jcvi.formats.fasta import get_first_rec, print_first_difference
from jcvi.apps.console import print_green
from jcvi.apps.base import ActionDispatcher, debug
debug()

myEmail = "htang@jcvi.org"
Entrez.email = myEmail


def batch_taxonomy(list_of_taxids):
    """
    Retrieve list of taxids, and generate latin names
    """
    for taxid in list_of_taxids:
        handle = Entrez.efetch(db='Taxonomy', id=taxid, retmode="xml")
        records = Entrez.read(handle)
        yield records[0]["ScientificName"]


def batch_entrez(list_of_terms, db="nucleotide", retmax=1, rettype="fasta"):
    """
    Retrieve multiple rather than a single record
    """

    for term in list_of_terms:

        logging.debug("search term %s" % term)
        success = False
        ids = None
        while not success:
            try:
                search_handle = Entrez.esearch(db=db, retmax=retmax, term=term)
                rec = Entrez.read(search_handle)
                success = True
                ids = rec["IdList"]
            except (urllib2.HTTPError, urllib2.URLError,
                    RuntimeError, KeyError) as e:
                logging.error(e)
                logging.debug("wait 5 seconds to reconnect...")
                time.sleep(5)

        assert ids

        for id in ids:
            success = False
            while not success:
                try:
                    fetch_handle = Entrez.efetch(db=db, id=id, rettype=rettype,
                            email=myEmail)
                    success = True
                except (urllib2.HTTPError, urllib2.URLError,
                        RuntimeError) as e:
                    logging.error(e)
                    logging.debug("wait 5 seconds to reconnect...")
                    time.sleep(5)

            yield id, term, fetch_handle


def main():

    actions = (
        ('fetch', 'fetch records from a list of GenBank accessions'),
        ('bisect', 'determine the version of the accession'),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bisect(args):
    """
    %prog bisect acc accession.fasta

    determine the version of the accession, based on a fasta file.
    This proceeds by a sequential search from xxxx.1 to the latest record.
    """
    p = OptionParser(bisect.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    acc, fastafile = args
    arec = get_first_rec(fastafile)

    valid = None
    for i in range(1, 100):
        term = "%s.%d" % (acc, i)
        try:
            query = list(batch_entrez([term]))
        except AssertionError as e:
            logging.debug("no records found for %s. terminating." % term)
            return

        id, term, handle = query[0]
        brec = SeqIO.parse(handle, "fasta").next()

        match = print_first_difference(arec, brec, ignore_case=True,
                ignore_N=True, rc=True)
        if match:
            valid = term
            break

    if valid:
        print
        print_green("%s matches the sequence in `%s`" % (valid, fastafile))


def fetch(args):
    """
    %prog fetch filename

    filename contains a list of terms to search
    """
    p = OptionParser(fetch.__doc__)

    valid_formats = ("fasta", "gb")
    p.add_option("--noversion", dest="noversion",
            default=False, action="store_true",
            help="Remove trailing accession versions")
    p.add_option("--format", default="fasta", choices=valid_formats,
            help="download format [default: %default]")
    p.add_option("--outdir", default=None,
            help="output directory, with accession number as filename")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    filename, = args
    if op.exists(filename):
        list_of_terms = [row.strip() for row in open(filename)]
        if opts.noversion:
            list_of_terms = [x.rsplit(".", 1)[0] for x in list_of_terms]
    else:
        # the filename is the search term
        list_of_terms = [filename.strip()]

    format = opts.format

    outfile = "{0}.{1}".format(filename.rsplit(".", 1)[0], format)
    if op.exists(outfile):
        logging.error("`{0}` found, overwrite (Y/N)?".format(outfile))
        while True:
            yesno = raw_input()
            if yesno == 'N':
                return
            elif yesno == 'Y':
                break
            else:
                logging.error("You have to answer Y/N ..")
                continue

    outdir = opts.outdir
    if outdir and not op.exists(outdir):
        logging.debug("`%s` not found, creating new." % outdir)
        os.mkdir(outdir)

    if not outdir:
        fw = open(outfile, "w")

    seen = set()
    for id, term, handle in batch_entrez(list_of_terms, rettype=format):
        rec = handle.read()
        if id in seen:
            logging.error("duplicate key (%s) found" % rec)
            continue

        if outdir:
            fw = open(op.join(outdir, term), "w")

        print >> fw, rec
        print >> fw

        seen.add(id)


if __name__ == '__main__':
    main()
