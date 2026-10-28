"""
Wrapper for fetching data from various online repositories \
(Entrez, Ensembl, Phytozome, and SRA)
"""

import json
import os
import os.path as op
from os.path import join as urljoin
import re
import shutil
import sys
import time
from typing import Dict, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from Bio import Entrez, SeqIO
from more_itertools import grouper

from ..formats.base import must_open
from ..formats.fasta import print_first_difference
from ..formats.fastq import fromsra
from ..utils.cbook import tile
from rich.console import Console

console = Console()

from .base import (
    ActionDispatcher,
    OptionParser,
    cleanup,
    download,
    get_email_address,
    logger,
    ls_ftp,
    mkdir,
)

myEmail = get_email_address()
Entrez.email = myEmail
JGI_FILE_LIST_URL = "https://files.jgi.doe.gov/phytozome_file_list/"
JGI_DOWNLOAD_URL = "https://files-download.jgi.doe.gov/download_files/{}/"
JGI_PORTAL_URL = "https://data.jgi.doe.gov/refine-download/phytozome"
JGI_TOKEN_ENV = "JGI_SESSION_TOKEN"
# Files on tape, need a restore request before download
JGI_ARCHIVED_STATUS = {"PURGED", "RESTORE_REGISTERED", "RESTORE_IN_PROGRESS"}


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

        logger.debug("Search term %s", term)
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
                logger.error(e)
                logger.debug("wait 5 seconds to reconnect...")
                time.sleep(5)

        if not ids:
            logger.error("term {0} not found".format(term))
            continue

        assert ids
        nids = len(ids)
        if nids > 1:
            logger.debug("A total of {0} results found.".format(nids))

        if batchsize != 1:
            logger.debug("Use a batch size of {0}.".format(batchsize))

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
                    logger.error(e)
                    logger.debug("wait 5 seconds to reconnect...")
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
    p.add_argument("--version", default="75", help="Ensembl version")
    opts, args = p.parse_args(args)

    version = opts.version
    url = "ftp://ftp.ensembl.org/pub/release-{0}/".format(version)
    fasta_url = url + "fasta/"

    valid_species = [x for x in ls_ftp(fasta_url) if "." not in x]
    doc = "\n".join((ensembl.__doc__, tile(valid_species)))
    p.usage = doc.replace("%prog", "%(prog)s")

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


def get_jgi_token() -> Optional[str]:
    """Read JGI session token from the environment, or prompt for it."""
    token = os.environ.get(JGI_TOKEN_ENV)
    if not token and console.is_terminal:
        token = console.input("[bold green]JGI session token: ", password=True)
    if not token:
        return None
    token = token.strip()
    return token if token.startswith("Bearer ") else f"Bearer {token}"


def query_phytozome(
    species: str, version: Optional[str] = None, assembly: bool = False
) -> List[dict]:
    """Query the JGI file list for a species name or a genome id."""
    pattern = r".*[._](gene\.gff3|cds\.fa)\.gz"
    if assembly:
        pattern += r"|.*\.fa\.gz"
    params = {"api_version": 2, "x": 50, "ff[file_name_pattern]": pattern}
    if species.isdigit():
        params["genome_id"] = species
    else:
        params["organism"] = species
        if version:
            params["version"] = version
    with urlopen(f"{JGI_FILE_LIST_URL}?{urlencode(params)}") as fp:
        return json.load(fp).get("organisms", [])


def classify_phytozome_file(f: dict) -> Optional[str]:
    """Label a file record as gff, cds or asm (unmasked assembly)."""
    name, meta = f["file_name"], f.get("metadata", {})
    if "lowConfidence" in name:
        return None
    if re.search(r"[._]gene\.gff3\.gz$", name):
        return "gff"
    if re.search(r"[._]cds\.fa\.gz$", name):
        return "cds"
    # metadata.content is unreliable (e.g. "unmashed"), check name instead
    if (
        name.endswith(".fa.gz")
        and meta.get("type") == "assembly"
        and "masked" not in name
    ):
        return "asm"
    return None


def _genome_id(organism: dict) -> int:
    return int(organism["id"].rsplit("-", 1)[-1])


def _rank_phytozome_file(f: dict, version: str) -> Tuple[bool, bool, str]:
    """Prefer files on disk, then named after the genome version, then newest."""
    return (
        f.get("file_status") not in JGI_ARCHIVED_STATUS,
        f"_{version}." in f["file_name"],
        f.get("added_date", ""),
    )


def select_phytozome_files(
    organisms: List[dict], assembly: bool = False
) -> Dict[str, dict]:
    """Pick files from the newest genome that has all requested file kinds."""
    kinds = ("gff", "cds", "asm") if assembly else ("gff", "cds")
    for o in sorted(organisms, key=_genome_id, reverse=True):
        candidates: Dict[Optional[str], List[dict]] = {}
        for f in o.get("files", []):
            candidates.setdefault(classify_phytozome_file(f), []).append(f)
        if not all(k in candidates for k in kinds):
            continue
        logger.info("Selected %s (%s)", o["name"], o["id"])
        if o.get("data_utilization_status") == "Restricted":
            logger.warning("%s is under restricted data use policy", o["name"])
        meta = o.get("top_hit", {}).get("metadata", {}).get("phytozome", {})
        version = meta.get("version", "")
        return {
            k: max(candidates[k], key=lambda f: _rank_phytozome_file(f, version))
            for k in kinds
        }
    return {}


def download_jgi_file(f: dict, token: str) -> Optional[str]:
    """Download a single file from JGI Data Portal, return the local path."""
    filename = f["file_name"]
    if op.exists(filename):
        logger.info("File `%s` exists. Download skipped.", filename)
        return filename
    req = Request(JGI_DOWNLOAD_URL.format(f["_id"]), headers={"Authorization": token})
    logger.info("Downloading `%s` (%s)", filename, f.get("file_status"))
    try:
        with urlopen(req) as src, open(filename, "wb") as dst:
            shutil.copyfileobj(src, dst)
    except HTTPError as e:
        cleanup(filename)
        if e.code in (401, 403):
            logger.error("Invalid or expired JGI session token")
        elif e.code == 409:
            logger.error(
                "`%s` is archived on tape, request a restore at %s",
                filename,
                JGI_PORTAL_URL,
            )
        else:
            logger.error("Failed to download `%s`: %s", filename, e)
        return None
    except URLError as e:
        cleanup(filename)
        logger.error("Failed to download `%s`: %s", filename, e)
        return None
    return filename


def phytozome(args):
    """
    %prog phytozome species

    Retrieve genomes and annotations from Phytozome using JGI Data Portal API.
    Use comma to give a list of species, which are Phytozome short names
    (the newest genome is picked) or genome ids (to pick a specific genome).
    Browse available genomes at https://phytozome-next.jgi.doe.gov/. For example:

    $ %prog phytozome Athaliana,Vvinifera,Osativa,Sbicolor,Slycopersicum

    Downloading requires a JGI session token. Log in at https://data.jgi.doe.gov,
    click your avatar -> "Copy My Session Token", then set it to the environment
    variable JGI_SESSION_TOKEN, or paste it when prompted.
    """
    p = OptionParser(phytozome.__doc__)
    p.add_argument(
        "--version",
        default="current",
        help="Phytozome release, e.g. 13, ignored for genome ids",
    )
    p.add_argument(
        "--assembly",
        default=False,
        action="store_true",
        help="Download assembly",
    )
    p.add_argument(
        "--format",
        default=False,
        action="store_true",
        help="Format to CDS and BED for synteny inference",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (species,) = args
    token = get_jgi_token()
    if token is None:
        logger.error("No JGI session token found, set %s", JGI_TOKEN_ENV)
        sys.exit(1)

    for s in species.split(","):
        res = download_species_phytozome(
            s, token, version=opts.version, assembly=opts.assembly
        )
        gff, cds = res.get("gff"), res.get("cds")
        if not (gff and cds):
            logger.error("No annotation downloaded for %s", s)
            continue
        if opts.format:
            format_bed_and_cds(s, gff, cds)


def download_species_phytozome(
    species: str,
    token: str,
    version: Optional[str] = None,
    assembly: bool = False,
) -> Dict[str, str]:
    """Download annotation GFF, CDS FASTA and optionally assembly FASTA.

    Returns a dict that maps gff|cds|asm to local file paths.
    """
    organisms = query_phytozome(species, version=version, assembly=assembly)
    files = select_phytozome_files(organisms, assembly=assembly)
    if not files:
        logger.error("%s not found in Phytozome (version=%s)", species, version)
        return {}
    res = {}
    for kind, f in files.items():
        path = download_jgi_file(f, token)
        if path:
            res[kind] = path
    return res


def phytozome9(args):
    """
    %prog phytozome9 species

    Retrieve genomes and annotations from phytozome FTP. Available species
    listed below. Use comma to give a list of species to download. For example:

    $ %prog phytozome9 Athaliana,Vvinifera,Osativa,Sbicolor,Slycopersicum
    """
    p = OptionParser(phytozome9.__doc__)
    p.add_argument(
        "--assembly",
        default=False,
        action="store_true",
        help="Download assembly",
    )
    p.add_argument(
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
    p.usage = doc.replace("%prog", "%(prog)s")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (species,) = args
    if species == "all":
        species = ",".join(valid_species)

    species = species.split(",")

    for s in species:
        res = download_species_phytozome9(s, valid_species, url, assembly=opts.assembly)
        if not res:
            logger.error("No files downloaded")
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
    from jcvi.formats.fasta import format as fasta_format
    from jcvi.formats.gff import bed as gff_bed

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
    logger.debug("Found magic number for {0}: {1}".format(species, magic))

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
        logger.debug(
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
            logger.debug(f"no records found for {term}. terminating. {e}")
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
        console.print()
        console.print("[green]{} matches the sequence in `{}`".format(valid, fastafile))


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

    p.add_argument(
        "--noversion",
        dest="noversion",
        default=False,
        action="store_true",
        help="Remove trailing accession versions",
    )
    p.add_argument(
        "--format",
        default="fasta",
        choices=valid_formats,
        help="download format",
    )
    p.add_argument(
        "--database",
        default="nuccore",
        choices=valid_databases,
        help="search database",
    )
    p.add_argument(
        "--retmax",
        default=1000000,
        type=int,
        help="how many results to return",
    )
    p.add_argument(
        "--skipcheck",
        default=False,
        action="store_true",
        help="turn off prompt to check file existence",
    )
    p.add_argument(
        "--batchsize",
        default=500,
        type=int,
        help="download the results in batch for speed-up",
    )
    p.set_outdir(outdir=None)
    p.add_argument("--outprefix", default="out", help="output file name prefix")
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
            logger.error("Duplicate key ({0}) found".format(rec))
            continue

        totalsize += size
        print(rec, file=fw)
        print(file=fw)

        seen.add(id)

    if seen:
        console.print(
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

    p.add_argument(
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
        logger.error(
            "Incorrect SRA identifier format "
            + "[should be like SRR126150, SRR1001901. "
            + "len(identifier) should be between 9-10 characters]"
        )
        sys.exit()

    prefix, subprefix = m.group(1), "{0}{1}".format(m.group(1), m.group(2))
    download_url = urljoin(
        sra_base_url, prefix, subprefix, term, "{0}.sra".format(term)
    )

    logger.debug("Downloading file: {0}".format(download_url))
    return download(download_url)


if __name__ == "__main__":
    main()
