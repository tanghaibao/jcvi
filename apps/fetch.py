"""
Wrapper for fetching data from various online repositories \
(Entrez, Ensembl, Phytozome, and SRA)
"""

import os.path as op
import sys
import time
import logging
import urllib2
import re
from os.path import join as urljoin

from Bio import Entrez, SeqIO

from jcvi.formats.base import must_open
from jcvi.formats.fasta import print_first_difference
from jcvi.utils.cbook import tile
from jcvi.utils.iter import grouper
from jcvi.apps.console import green
from jcvi.apps.base import OptionParser, ActionDispatcher, get_email_address, \
            mkdir, ls_ftp, download


myEmail = get_email_address()
Entrez.email = myEmail


def batch_taxonomy(list_of_taxids):
    """
    Convert list of taxids to Latin names
    """
    for taxid in list_of_taxids:
        handle = Entrez.efetch(db='Taxonomy', id=taxid, retmode="xml")
        records = Entrez.read(handle)
        yield records[0]["ScientificName"]


def batch_taxids(list_of_names):
    """
    Opposite of batch_taxonomy():

    Convert list of Latin names to taxids
    """
    for name in list_of_names:
        handle = Entrez.esearch(db='Taxonomy', term=name, retmode="xml")
        records = Entrez.read(handle)
        yield records["IdList"][0]


def batch_entrez(list_of_terms, db="nuccore", retmax=1, rettype="fasta",
            batchsize=1, email=myEmail):
    """
    Retrieve multiple rather than a single record
    """

    for term in list_of_terms:

        logging.debug("Search term %s" % term)
        success = False
        ids = None
        if not term:
            continue

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

        if not ids:
            logging.error("term {0} not found".format(term))
            continue

        assert ids
        nids = len(ids)
        if nids > 1:
            logging.debug("A total of {0} results found.".format(nids))

        if batchsize != 1:
            logging.debug("Use a batch size of {0}.".format(batchsize))

        ids = list(grouper(ids, batchsize))

        for id in ids:
            id = [x for x in id if x]
            size = len(id)
            id = ",".join(id)

            success = False
            while not success:
                try:
                    fetch_handle = Entrez.efetch(db=db, id=id, rettype=rettype,
                            email=email)
                    success = True
                except (urllib2.HTTPError, urllib2.URLError,
                        RuntimeError) as e:
                    logging.error(e)
                    logging.debug("wait 5 seconds to reconnect...")
                    time.sleep(5)

            yield id, size, term, fetch_handle


def main():

    actions = (
        ('entrez', 'fetch records from entrez using a list of GenBank accessions'),
        ('bisect', 'determine the version of the accession by querying entrez'),
        ('phytozome', 'retrieve genomes and annotations from phytozome'),
        ('ensembl', 'retrieve genomes and annotations from ensembl'),
        ('sra', 'retrieve files from SRA via the sra-instant FTP'),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def ensembl(args):
    """
    %prog ensembl species

    Retrieve genomes and annotations from ensembl FTP. Available species
    listed below. Use comma to give a list of species to download. For example:

    $ %prog ensembl danio_rerio,gasterosteus_aculeatus
    """
    p = OptionParser(ensembl.__doc__)
    p.add_option("--version", default="75",
                 help="Ensembl version [default: %default]")
    opts, args = p.parse_args(args)

    version = opts.version
    url = "ftp://ftp.ensembl.org/pub/release-{0}/".format(version)
    fasta_url = url + "fasta/"

    valid_species = [x for x in ls_ftp(fasta_url) if "." not in x]
    doc = "\n".join((ensembl.__doc__, tile(valid_species)))
    p.set_usage(doc)

    if len(args) != 1:
        sys.exit(not p.print_help())

    species, = args
    species = species.split(",")
    for s in species:
        download_species_ensembl(s, valid_species, url)


def download_species_ensembl(species, valid_species, url):
    assert species in valid_species, \
            "{0} is not in the species list".format(species)

    # We want to download assembly and annotation for given species
    ann_url = urljoin(url, "gtf/{0}".format(species))
    cds_url = urljoin(url, "fasta/{0}/cds".format(species))

    for u in (ann_url, cds_url):
        valid_files = [x for x in ls_ftp(u) if x.endswith(".gz")]
        for f in valid_files:
            f = urljoin(u, f)
            download(f)


def phytozome(args):
    """
    %prog phytozome species

    Retrieve genomes and annotations from phytozome FTP. Available species
    listed below. Use comma to give a list of species to download. For example:

    $ %prog phytozome Athaliana,Vvinifera,Osativa,Sbicolor,Slycopersicum
    """
    p = OptionParser(phytozome.__doc__)
    p.add_option("--version", default="9.0",
                 help="Phytozome version [default: %default]")
    p.add_option("--assembly", default=False, action="store_true",
                 help="Download assembly [default: %default]")
    opts, args = p.parse_args(args)

    url = "ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v{0}/".\
                    format(opts.version)
    valid_species = [x for x in ls_ftp(url) if "." not in x]

    doc = "\n".join((phytozome.__doc__, tile(valid_species)))
    p.set_usage(doc)

    if len(args) != 1:
        sys.exit(not p.print_help())

    species, = args
    species = species.split(",")
    for s in species:
        download_species_phytozome(s, valid_species, url, assembly=opts.assembly)


def download_species_phytozome(species, valid_species, url, assembly=False):
    from os.path import join as urljoin

    assert species in valid_species, \
            "{0} is not in the species list".format(species)

    # We want to download assembly and annotation for given species
    surl = urljoin(url, species)
    contents = [x for x in ls_ftp(surl) if x.endswith("_readme.txt")]
    magic = contents[0].split("_")[1]  # Get the magic number
    logging.debug("Found magic number for {0}: {1}".format(species, magic))

    pf = "{0}_{1}".format(species, magic)
    asm_url = urljoin(surl, "assembly/{0}.fa.gz".format(pf))
    ann_url = urljoin(surl, "annotation/{0}_gene.gff3.gz".format(pf))
    cds_url = urljoin(surl, "annotation/{0}_cds.fa.gz".format(pf))
    if assembly:
        download(asm_url)
    for u in (ann_url, cds_url):
        download(u)


def get_first_rec(fastafile):
    """
    Returns the first record in the fastafile
    """
    f = list(SeqIO.parse(fastafile, "fasta"))

    if len(f) > 1:
        logging.debug("{0} records found in {1}, using the first one".\
                format(len(f), fastafile))

    return f[0]


def bisect(args):
    """
    %prog bisect acc accession.fasta

    determine the version of the accession by querying entrez, based on a fasta file.
    This proceeds by a sequential search from xxxx.1 to the latest record.
    """
    p = OptionParser(bisect.__doc__)
    p.set_email()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    acc, fastafile = args
    arec = get_first_rec(fastafile)

    valid = None
    for i in range(1, 100):
        term = "%s.%d" % (acc, i)
        try:
            query = list(batch_entrez([term], email=opts.email))
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
        print green("%s matches the sequence in `%s`" % (valid, fastafile))


def entrez(args):
    """
    %prog entrez <filename|term>

    `filename` contains a list of terms to search. Or just one term. If the
    results are small in size, e.g. "--format=acc", use "--batchsize=100" to speed
    the download.
    """
    p = OptionParser(entrez.__doc__)

    allowed_databases = {"fasta": ["genome", "nuccore", "nucgss", "protein", "nucest"],
                         "asn.1": ["genome", "nuccore", "nucgss", "protein"],
                         "gb"   : ["genome", "nuccore", "nucgss"],
                         "est"  : ["nucest"],
                         "gss"  : ["nucgss"],
                         "acc"  : ["nuccore"],
                        }

    valid_formats = tuple(allowed_databases.keys())
    valid_databases = ("genome", "nuccore", "nucest", "nucgss", "protein")

    p.add_option("--noversion", dest="noversion",
            default=False, action="store_true",
            help="Remove trailing accession versions")
    p.add_option("--format", default="fasta", choices=valid_formats,
            help="download format [default: %default]")
    p.add_option("--database", default="nuccore", choices=valid_databases,
            help="search database [default: %default]")
    p.add_option("--retmax", default=1000000, type="int",
            help="how many results to return [default: %default]")
    p.add_option("--skipcheck", default=False, action="store_true",
            help="turn off prompt to check file existence [default: %default]")
    p.add_option("--batchsize", default=500, type="int",
            help="download the results in batch for speed-up [default: %default]")
    p.add_option("--outdir", default=None,
            help="output directory, with accession number as filename")
    p.add_option("--outprefix", default="out",
            help="output file name prefix [default: %default]")
    p.set_email()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    filename, = args
    if op.exists(filename):
        pf = filename.rsplit(".", 1)[0]
        list_of_terms = [row.strip() for row in open(filename)]
        if opts.noversion:
            list_of_terms = [x.rsplit(".", 1)[0] for x in list_of_terms]
    else:
        pf = filename
        # the filename is the search term
        list_of_terms = [filename.strip()]

    fmt = opts.format
    database = opts.database
    batchsize = opts.batchsize

    assert database in allowed_databases[fmt], \
        "For output format '{0}', allowed databases are: {1}".\
        format(fmt, allowed_databases[fmt])
    assert batchsize >= 1, "batchsize must >= 1"

    if " " in pf:
        pf = opts.outprefix

    outfile = "{0}.{1}".format(pf, fmt)

    outdir = opts.outdir
    if outdir:
        mkdir(outdir)

    # If noprompt, will not check file existence
    if not outdir:
        fw = must_open(outfile, "w", checkexists=True, \
                skipcheck=opts.skipcheck)
        if fw is None:
            return

    seen = set()
    totalsize = 0
    for id, size, term, handle in batch_entrez(list_of_terms, retmax=opts.retmax, \
                                 rettype=fmt, db=database, batchsize=batchsize, \
                                 email=opts.email):
        if outdir:
            outfile = urljoin(outdir, "{0}.{1}".format(term, fmt))
            fw = must_open(outfile, "w", checkexists=True, \
                    skipcheck=opts.skipcheck)
            if fw is None:
                continue

        rec = handle.read()
        if id in seen:
            logging.error("Duplicate key ({0}) found".format(rec))
            continue

        totalsize += size
        print >> fw, rec
        print >> fw

        seen.add(id)

    if seen:
        print >> sys.stderr, "A total of {0} {1} records downloaded.".\
                format(totalsize, fmt.upper())

    return outfile


def sra(args):
    """
    %prog sra term

    Given an SRA run ID, fetch the corresponding .sra file
    from the sra-instant FTP
    """
    sra_base_url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/"
    sra_run_id_re = re.compile(r'^([DES]{1}RR)(\d{3})(\d{3,4})$')

    p = OptionParser(sra.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    term, = args

    m = re.search(sra_run_id_re, term)
    if m is None:
        logging.error("Incorrect SRA identifier format " + \
                "[should be like SRR126150, SRR1001901. " + \
                "len(identifier) should be between 9-10 characters]")
        sys.exit()

    prefix, subprefix = m.group(1), "{0}{1}".format(m.group(1), m.group(2))
    download_url = urljoin(sra_base_url, prefix, subprefix, term, "{0}.sra".format(term))

    logging.debug("Downloading file: {0}".format(download_url))
    download(download_url)


if __name__ == '__main__':
    main()
