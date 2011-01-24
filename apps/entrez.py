"""
Wrapper for calling Bio.Entrez tools to get the sequence from a list of IDs
"""

import sys
import time
import logging
import urllib2

from optparse import OptionParser
from Bio import Entrez

from jcvi.apps.base import debug


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
        search_handle = Entrez.esearch(db=db, retmax=retmax, term=term)
        rec = Entrez.read(search_handle)
        ids = rec["IdList"]

        if not ids:
            logging.error("term %s not found in db %s" % (term, db))

        for id in ids:
            success = False 
            while not success:
                try:
                    fetch_handle = Entrez.efetch(db=db, id=id, rettype=rettype,
                            email=myEmail)
                    success = True
                except (urllib2.HTTPError, urllib2.URLError) as e:
                    logging.error(str(e))
                    logging.debug("wait 5 seconds to reconnect...")
                    time.sleep(5)

            yield id, fetch_handle.read()


def main():
    """
    %prog filename

    filename contains a list of terms to search 
    """
    debug()

    p = OptionParser(main.__doc__)

    valid_formats = ("fasta", "gb")
    p.add_option("--format", default="fasta", choices=valid_formats,
            help="download format [default: %default]")
    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    filename = args[0]
    list_of_terms = [row.strip() for row in open(filename)]
    
    seen = set()
    for id, rec in batch_entrez(list_of_terms, rettype=opts.format):
        if id in seen:
            logging.error("duplicate key (%s) found" % rec)
            continue
        else:
            print rec

        seen.add(id)


if __name__ == '__main__':
    main()
