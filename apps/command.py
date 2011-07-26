"""
Commonly performed commands.
"""

import shutil
import logging

from jcvi.utils.cbook import depends
from jcvi.apps.base import sh


# Install locations of common binaries
BDPATH = "~/bin/"
CAPATH = "~/bin/Linux-amd64/bin/"
CDPATH="~/htang/export/cd-hit-v4.5.5-2011-03-31/"


@depends
def run_formatdb(infile=None, outfile=None):
    cmd = "formatdb -i {0} -p F".format(infile)
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

    cmd = 'blastall -p blastn -i {0}'.format(infile)
    cmd += ' -d {0} -q -5 -G 3 -E 3 -F "m D"'.format(db)
    cmd += ' -e 0.01 -Y 1.75e12 -m 8 -o {0} -a 8'.format(outfile)
    sh(cmd)


@depends
def run_megablast(infile=None, outfile=None, db=None, pctid=98, hitlen=100):
    nin = db + ".nin"
    run_formatdb(infile=db, outfile=nin)

    cmd = "megablast -i {0} -d {1}".format(infile, db)
    cmd += " -e 0.001 -m 8 -o {0} -a 8".format(outfile)
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
