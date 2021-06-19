"""
Wrapper for fetching data from various online repositories \
(Entrez, Ensembl, Phytozome, and SRA)
"""
import logging
import os.path as op
import re
import sys
import time

from os.path import join as urljoin
from urllib.error import HTTPError, URLError

from Bio import Entrez, SeqIO
from more_itertools import grouper

from jcvi.formats.base import FileShredder, must_open
from jcvi.formats.fasta import print_first_difference
from jcvi.formats.fastq import fromsra
from jcvi.utils.cbook import tile
from jcvi.utils.console import printf
from jcvi.apps.base import (
    OptionParser,
    ActionDispatcher,
    get_email_address,
    mkdir,
    ls_ftp,
    download,
    sh,
    last_updated,
    which,
)


myEmail = get_email_address()
Entrez.email = myEmail
PHYTOZOME_COOKIES = ".phytozome_cookies"


def batch_taxonomy(list_of_taxids):
    """
    Convert list of taxids to Latin names
    """
    for taxid in list_of_taxids:
        handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        yield records[0]["ScientificName"]


def batch_taxids(list_of_names):
    """
    Opposite of batch_taxonomy():

    Convert list of Latin names to taxids
    """
    for name in list_of_names:
        handle = Entrez.esearch(db="Taxonomy", term=name, retmode="xml")
        records = Entrez.read(handle)
        yield records["IdList"][0]


def batch_entrez(
    list_of_terms, db="nuccore", retmax=1, rettype="fasta", batchsize=1, email=myEmail
):
    """
    Retrieve multiple rather than a single record
    """

    for term in list_of_terms:

        logging.debug("Search term %s", term)
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
            except (HTTPError, URLError, RuntimeError, KeyError) as e:
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
                    fetch_handle = Entrez.efetch(
                        db=db, id=id, rettype=rettype, email=email
                    )
                    success = True
                except (HTTPError, URLError, RuntimeError) as e:
                    logging.error(e)
                    logging.debug("wait 5 seconds to reconnect...")
                    time.sleep(5)

            yield id, size, term, fetch_handle


def main():

    actions = (
        ("entrez", "fetch records from entrez using a list of GenBank accessions"),
        ("bisect", "determine the version of the accession by querying entrez"),
        (
            "phytozome9",
            "retrieve genomes and annotations from phytozome version 9.0 (legacy)",
        ),
        ("phytozome", "retrieve genomes and annotations from phytozome"),
        ("ensembl", "retrieve genomes and annotations from ensembl"),
        ("sra", "retrieve files from SRA via the sra-instant FTP"),
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
    p.add_option("--version", default="75", help="Ensembl version")
    opts, args = p.parse_args(args)

    version = opts.version
    url = "ftp://ftp.ensembl.org/pub/release-{0}/".format(version)
    fasta_url = url + "fasta/"

    valid_species = [x for x in ls_ftp(fasta_url) if "." not in x]
    doc = "\n".join((ensembl.__doc__, tile(valid_species)))
    p.set_usage(doc)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (species,) = args
    species = species.split(",")
    for s in species:
        download_species_ensembl(s, valid_species, url)


def download_species_ensembl(species, valid_species, url):
    assert species in valid_species, "{0} is not in the species list".format(species)

    # We want to download assembly and annotation for given species
    ann_url = urljoin(url, "gtf/{0}".format(species))
    cds_url = urljoin(url, "fasta/{0}/cds".format(species))

    for u in (ann_url, cds_url):
        valid_files = [x for x in ls_ftp(u) if x.endswith(".gz")]
        for f in valid_files:
            f = urljoin(u, f)
            download(f)


def get_cookies(cookies=PHYTOZOME_COOKIES):
    from jcvi.utils.console import console

    # Check if cookies is still good
    if op.exists(cookies) and last_updated(cookies) < 3600:
        return cookies

    if console.is_terminal:
        username = console.input("[bold green]Phytozome Login: ")
        pw = console.input("[bold green]Phytozome Password: ", password=True)
    else:
        username, pw = None, None
    curlcmd = which("curl")
    if curlcmd is None:
        logging.error("curl command not installed. Aborting.")
        return None
    cmd = "{} https://signon.jgi.doe.gov/signon/create".format(curlcmd)
    cmd += " --data-urlencode 'login={0}' --data-urlencode 'password={1}' -b {2} -c {2}".format(
        username, pw, cookies
    )
    sh(cmd, outfile="/dev/null", errfile="/dev/null", log=False)
    if not op.exists(cookies):
        logging.error("Cookies file `{}` not created. Aborting.".format(cookies))
        return None

    return cookies


def phytozome(args):
    """
    %prog phytozome species

    Retrieve genomes and annotations from phytozome using Globus API. Available
    species listed below. Use comma to give a list of species to download. For
    example:

    $ %prog phytozome Athaliana,Vvinifera,Osativa,Sbicolor,Slycopersicum

    The downloader will prompt you to enter Phytozome user name and password
    during downloading. Please register for a login at:
    https://phytozome.jgi.doe.gov/pz/portal.html.
    """
    from jcvi.apps.biomart import GlobusXMLParser

    p = OptionParser(phytozome.__doc__)
    p.add_option(
        "--version",
        default="12",
        choices=("9", "10", "11", "12", "12_unrestricted", "13"),
        help="Phytozome version",
    )
    p.add_option(
        "--assembly",
        default=False,
        action="store_true",
        help="Download assembly",
    )
    p.add_option(
        "--format",
        default=False,
        action="store_true",
        help="Format to CDS and BED for synteny inference",
    )
    p.set_downloader()
    opts, args = p.parse_args(args)

    downloader = opts.downloader
    directory_listing = ".phytozome_directory_V{}.xml".format(opts.version)
    # Get directory listing
    base_url = "http://genome.jgi.doe.gov"
    dlist = "{}/ext-api/downloads/get-directory?organism=PhytozomeV{}".format(
        base_url, opts.version
    )

    # Make sure we have a valid cookies
    cookies = get_cookies()
    if cookies is None:
        logging.error("Error fetching cookies ... cleaning up")
        FileShredder([directory_listing])
        sys.exit(1)

    # Proceed to use the cookies and download the species list
    try:
        download(
            dlist,
            filename=directory_listing,
            cookies=cookies,
            downloader=downloader,
        )
        g = GlobusXMLParser(directory_listing)
    except:
        logging.error("Error downloading directory listing ... cleaning up")
        FileShredder([directory_listing, cookies])
        sys.exit(1)

    genomes = g.get_genomes()
    valid_species = genomes.keys()
    species_tile = tile(valid_species)
    p.set_usage("\n".join((phytozome.__doc__, species_tile)))

    if len(args) != 1:
        sys.exit(not p.print_help())

    (species,) = args
    if species == "all":
        species = ",".join(valid_species)

    species = species.split(",")
    for s in species:
        res = download_species_phytozome(
            genomes,
            s,
            valid_species,
            base_url,
            cookies,
            assembly=opts.assembly,
            downloader=downloader,
        )
        if not res:
            logging.error("No files downloaded")
        gff, fa = res.get("gff"), res.get("cds")
        if opts.format:
            format_bed_and_cds(s, gff, fa)


def download_species_phytozome(
    genomes, species, valid_species, base_url, cookies, assembly=False, downloader=None
):
    """Download assembly FASTA and annotation GFF.

    Args:
        genomes (dict): Dictionary parsed from Globus XML.
        species (str): Target species to download.
        valid_species (List[str]): Allowed set of species
        base_url (str): URL.
        cookies (str): cookies file path.
        assembly (bool, optional): Do we download assembly FASTA (can be big).
        Defaults to False.
        downloader (str, optional): Use a given downloader. One of wget|curl|powershell|insecure.
        Defaults to None.
    """
    assert species in valid_species, "{} is not in the species list".format(species)
    res = {}
    genome = genomes.get(species)
    if not genome:
        return res

    genome_assembly = genome.get("assembly")
    if assembly and genome_assembly:
        asm_name = next(x for x in genome_assembly if x.endswith(".fa.gz"))
        if asm_name:
            res["asm"] = genome_assembly.download(
                asm_name, base_url, cookies, downloader=downloader
            )

    genome_annotation = genome.get("annotation")
    if genome_annotation:
        gff_name = next(x for x in genome_annotation if x.endswith(".gene.gff3.gz"))
        if gff_name:
            res["gff"] = genome_annotation.download(
                gff_name, base_url, cookies, downloader=downloader
            )
        cds_name = next(x for x in genome_annotation if x.endswith(".cds.fa.gz"))
        if cds_name:
            res["cds"] = genome_annotation.download(
                cds_name, base_url, cookies, downloader=downloader
            )

    return res


def phytozome9(args):
    """
    %prog phytozome9 species

    Retrieve genomes and annotations from phytozome FTP. Available species
    listed below. Use comma to give a list of species to download. For example:

    $ %prog phytozome9 Athaliana,Vvinifera,Osativa,Sbicolor,Slycopersicum
    """
    p = OptionParser(phytozome9.__doc__)
    p.add_option(
        "--assembly",
        default=False,
        action="store_true",
        help="Download assembly",
    )
    p.add_option(
        "--format",
        default=False,
        action="store_true",
        help="Format to CDS and BED for synteny inference",
    )
    opts, args = p.parse_args(args)

    version = "9.0"
    url = "ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v{0}/".format(version)
    valid_species = [x for x in ls_ftp(url) if "." not in x]

    doc = "\n".join((phytozome9.__doc__, tile(valid_species)))
    p.set_usage(doc)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (species,) = args
    if species == "all":
        species = ",".join(valid_species)

    species = species.split(",")

    for s in species:
        res = download_species_phytozome9(s, valid_species, url, assembly=opts.assembly)
        if not res:
            logging.error("No files downloaded")
        gff, cdsfa = res.get("gff"), res.get("cds")
        if opts.format:
            format_bed_and_cds(s, gff, cdsfa)


def format_bed_and_cds(species, gff, cdsfa):
    """Run gff.format() and fasta.format() to generate BED and CDS files.
    This prepares the input files for the MCscan synteny workflow.

    https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)

    Args:
        species (str): Name of the species
        gff (str): Path to the GFF file
        cdsfa (str): Path to the FASTA file
    """
    from jcvi.formats.gff import bed as gff_bed
    from jcvi.formats.fasta import format as fasta_format

    # We have to watch out when the gene names and mRNA names mismatch, in which
    # case we just extract the mRNA names
    use_IDs = set()
    use_mRNAs = {
        "Cclementina",
        "Creinhardtii",
        "Csinensis",
        "Fvesca",
        "Lusitatissimum",
        "Mesculenta",
        "Mguttatus",
        "Ppersica",
        "Pvirgatum",
        "Rcommunis",
        "Sitalica",
        "Tcacao",
        "Thalophila",
        "Vcarteri",
        "Vvinifera",
        "Zmays",
    }
    key = "ID" if species in use_IDs else "Name"
    ttype = "mRNA" if species in use_mRNAs else "gene"
    bedfile = species + ".bed"
    cdsfile = species + ".cds"
    gff_bed([gff, "--type={}".format(ttype), "--key={}".format(key), "-o", bedfile])
    fasta_format([cdsfa, cdsfile, r"--sep=|"])


def download_species_phytozome9(species, valid_species, base_url, assembly=False):
    assert species in valid_species, "{} is not in the species list".format(species)

    # We want to download assembly and annotation for given species
    surl = urljoin(base_url, species)
    contents = [x for x in ls_ftp(surl) if x.endswith("_readme.txt")]
    magic = contents[0].split("_")[1]  # Get the magic number
    logging.debug("Found magic number for {0}: {1}".format(species, magic))

    pf = "{0}_{1}".format(species, magic)
    asm_url = urljoin(surl, "assembly/{0}.fa.gz".format(pf))
    ann_url = urljoin(surl, "annotation/{0}_gene.gff3.gz".format(pf))
    cds_url = urljoin(surl, "annotation/{0}_cds.fa.gz".format(pf))
    res = {}
    if assembly:
        res["asm"] = download(asm_url)
    res["gff"] = download(ann_url)
    res["cds"] = download(cds_url)
    return res


def get_first_rec(fastafile):
    """
    Returns the first record in the fastafile
    """
    f = list(SeqIO.parse(fastafile, "fasta"))

    if len(f) > 1:
        logging.debug(
            "{0} records found in {1}, using the first one".format(len(f), fastafile)
        )

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
            logging.debug(f"no records found for {term}. terminating. {e}")
            return

        id, term, handle = query[0]
        brec = next(SeqIO.parse(handle, "fasta"))

        match = print_first_difference(
            arec, brec, ignore_case=True, ignore_N=True, rc=True
        )
        if match:
            valid = term
            break

    if valid:
        printf()
        printf("[green]{} matches the sequence in `{}`".format(valid, fastafile))


def entrez(args):
    """
    %prog entrez <filename|term>

    `filename` contains a list of terms to search. Or just one term. If the
    results are small in size, e.g. "--format=acc", use "--batchsize=100" to speed
    the download.
    """
    p = OptionParser(entrez.__doc__)

    allowed_databases = {
        "fasta": ["genome", "nuccore", "nucgss", "protein", "nucest"],
        "asn.1": ["genome", "nuccore", "nucgss", "protein", "gene"],
        "xml": ["genome", "nuccore", "nucgss", "nucest", "gene"],
        "gb": ["genome", "nuccore", "nucgss"],
        "est": ["nucest"],
        "gss": ["nucgss"],
        "acc": ["nuccore"],
    }

    valid_formats = tuple(allowed_databases.keys())
    valid_databases = ("genome", "nuccore", "nucest", "nucgss", "protein", "gene")

    p.add_option(
        "--noversion",
        dest="noversion",
        default=False,
        action="store_true",
        help="Remove trailing accession versions",
    )
    p.add_option(
        "--format",
        default="fasta",
        choices=valid_formats,
        help="download format",
    )
    p.add_option(
        "--database",
        default="nuccore",
        choices=valid_databases,
        help="search database",
    )
    p.add_option(
        "--retmax",
        default=1000000,
        type="int",
        help="how many results to return",
    )
    p.add_option(
        "--skipcheck",
        default=False,
        action="store_true",
        help="turn off prompt to check file existence",
    )
    p.add_option(
        "--batchsize",
        default=500,
        type="int",
        help="download the results in batch for speed-up",
    )
    p.set_outdir(outdir=None)
    p.add_option("--outprefix", default="out", help="output file name prefix")
    p.set_email()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (filename,) = args
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

    assert (
        database in allowed_databases[fmt]
    ), "For output format '{0}', allowed databases are: {1}".format(
        fmt, allowed_databases[fmt]
    )
    assert batchsize >= 1, "batchsize must >= 1"

    if " " in pf:
        pf = opts.outprefix

    outfile = "{0}.{1}".format(pf, fmt)

    outdir = opts.outdir
    if outdir:
        mkdir(outdir)

    # If noprompt, will not check file existence
    if not outdir:
        fw = must_open(outfile, "w", checkexists=True, skipcheck=opts.skipcheck)
        if fw is None:
            return

    seen = set()
    totalsize = 0
    for id, size, term, handle in batch_entrez(
        list_of_terms,
        retmax=opts.retmax,
        rettype=fmt,
        db=database,
        batchsize=batchsize,
        email=opts.email,
    ):
        if outdir:
            outfile = urljoin(outdir, "{0}.{1}".format(term, fmt))
            fw = must_open(outfile, "w", checkexists=True, skipcheck=opts.skipcheck)
            if fw is None:
                continue

        rec = handle.read()
        if id in seen:
            logging.error("Duplicate key ({0}) found".format(rec))
            continue

        totalsize += size
        print(rec, file=fw)
        print(file=fw)

        seen.add(id)

    if seen:
        printf(
            "A total of {0} {1} records downloaded.".format(totalsize, fmt.upper()),
        )

    return outfile


def sra(args):
    """
    %prog sra [term|term.ids]

    Given an SRA run ID, fetch the corresponding .sra file from the sra-instant FTP.
    The term can also be a file containing list of SRR ids, one per line.

    Once downloaded, the SRA file is processed through `fastq-dump` to produce
    FASTQ formatted sequence files, which are gzipped by default.
    """
    p = OptionParser(sra.__doc__)

    p.add_option(
        "--nogzip",
        dest="nogzip",
        default=False,
        action="store_true",
        help="Do not gzip the FASTQ generated by fastq-dump",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (term,) = args
    if op.isfile(term):
        terms = [x.strip() for x in open(term)]
    else:
        terms = [term]

    for term in terms:
        srafile = download_srr_term(term)
        pf = srafile.split(".")[0]
        mkdir(pf)
        _opts = [srafile, "--paired", "--outdir={0}".format(pf)]
        if not opts.nogzip:
            _opts.append("--compress=gzip")
        fromsra(_opts)


def download_srr_term(term):
    sra_base_url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/"
    sra_run_id_re = re.compile(r"^([DES]RR)(\d{3})(\d{3,4})$")

    m = re.search(sra_run_id_re, term)
    if m is None:
        logging.error(
            "Incorrect SRA identifier format "
            + "[should be like SRR126150, SRR1001901. "
            + "len(identifier) should be between 9-10 characters]"
        )
        sys.exit()

    prefix, subprefix = m.group(1), "{0}{1}".format(m.group(1), m.group(2))
    download_url = urljoin(
        sra_base_url, prefix, subprefix, term, "{0}.sra".format(term)
    )

    logging.debug("Downloading file: {0}".format(download_url))
    return download(download_url)


if __name__ == "__main__":
    main()
