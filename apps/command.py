"""
Commonly performed commands.
"""

import os
import shutil
import logging

from jcvi.utils.cbook import depends
from jcvi.apps.base import sh


# Install locations of common binaries
BDPATH = "~/bin/"
CAPATH = "~/bin/Linux-amd64/bin/"
CDPATH = "~/htang/export/cd-hit-v4.5.5-2011-03-31/"
BLPATH = "~/scratch/bin/"
JAVA = "/usr/local/bin/java-1.6.0 -Xmx4g"

@depends
def run_formatdb(infile=None, outfile=None):
    cmd = BLPATH + "makeblastdb -dbtype nucl -in {0}".format(infile)
    sh(cmd)


@depends
def run_blat(infile=None, outfile=None, db="UniVec_Core", pctid=95, hitlen=50):
    cmd = 'blat {0} {1} -out=blast8 {2}'.format(db, infile, outfile)
    sh(cmd)

    blatfile = outfile
    filtered_blatfile = outfile + ".P{0}L{1}".format(pctid, hitlen)
    run_blast_filter(infile=blatfile, outfile=filtered_blatfile,
            pctid=pctid, hitlen=hitlen)
    shutil.move(filtered_blatfile, blatfile)


@depends
def run_vecscreen(infile=None, outfile=None, db="UniVec_Core",
        pctid=None, hitlen=None):
    """
    BLASTN parameters reference:
    http://www.ncbi.nlm.nih.gov/VecScreen/VecScreen_docs.html
    """
    nin = db + ".nin"
    run_formatdb(infile=db, outfile=nin)

    cmd = BLPATH + "blastn -task blastn"
    cmd += " -query {0} -db {1} -out {2}".format(infile, db, outfile)
    cmd += " -penalty -5 -gapopen 4 -gapextend 4 -dust yes -soft_masking true"
    cmd += " -searchsp 1750000000000 -evalue 0.01 -outfmt 6 -num_threads 8"
    sh(cmd)


@depends
def run_megablast(infile=None, outfile=None, db=None, pctid=98, hitlen=100):
    nin = db + ".nin"
    run_formatdb(infile=db, outfile=nin)

    cmd = BLPATH + "blastn"
    cmd += " -query {0} -db {1} -out {2}".format(infile, db, outfile)
    cmd += " -evalue 0.01 -outfmt 6 -num_threads 8"
    sh(cmd)

    blastfile = outfile
    filtered_blastfile = outfile + ".P{0}L{1}".format(pctid, hitlen)
    run_blast_filter(infile=blastfile, outfile=filtered_blastfile,
            pctid=pctid, hitlen=hitlen)
    shutil.move(filtered_blastfile, blastfile)


def run_blast_filter(infile=None, outfile=None, pctid=95, hitlen=50):
    from jcvi.formats.blast import filter

    logging.debug("Filter BLAST result (pctid={0}, hitlen={1})".\
            format(pctid, hitlen))
    pctidopt = "--pctid={0}".format(pctid)
    hitlenopt = "--hitlen={0}".format(hitlen)
    filter([infile, pctidopt, hitlenopt])
