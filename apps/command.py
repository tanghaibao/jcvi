"""
Commonly performed commands.
"""

import logging

from jcvi.utils.cbook import depends


# Install locations of common binaries
BDPATH = "~/bin/"
CAPATH = "~/bin/Linux-amd64/bin/"
CDPATH="~/htang/export/cd-hit-v4.5.5-2011-03-31/"


@depends
def run_formatdb(infile=None, outfile=None):
    cmd = "formatdb -i {0} -p F".format(infile)
    sh(cmd)


def run_blast_filter(infile=None, outfile=None, pctid=95, hitlen=50):
    from jcvi.formats.blast import filter

    logging.debug("Filter BLAST result (pctid={0}, hitlen={1})".\
            format(pctid, hitlen))
    pctidopt = "--pctid={0}".format(pctid)
    hitlenopt = "--hitlen={0}".format(hitlen)
    filter([infile, pctidopt, hitlenopt])
