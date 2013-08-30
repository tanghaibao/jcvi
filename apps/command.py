"""
Commonly performed commands.
"""

import os.path as op
import sys
import shutil
import logging
import ConfigParser

from functools import partial

from jcvi.apps.base import debug, sh, is_exe, which
debug()


def getpath(cmd, name=None, url=None, cfg="~/.jcvirc", warn="exit"):
    """
    Get install locations of common binaries
    First, check ~/.jcvirc file to get the full path
    If not present, ask on the console and store
    """
    p = which(cmd)  # if in PATH, just returns it
    if p:
        return p

    PATH = "Path"
    config = ConfigParser.RawConfigParser()
    cfg = op.expanduser(cfg)
    changed = False
    if op.exists(cfg):
        config.read(cfg)

    assert name is not None, "Need a program name"

    try:
        fullpath = config.get(PATH, name)
    except ConfigParser.NoSectionError:
        config.add_section(PATH)
        changed = True
    except:
        pass

    try:
        fullpath = config.get(PATH, name)
    except ConfigParser.NoOptionError:
        msg = "=== Configure path for {0} ===\n".format(name, cfg)
        if url:
            msg += "URL: {0}\n".format(url)
        msg += "[Directory that contains `{0}`]: ".format(cmd)
        fullpath = raw_input(msg).strip()
        config.set(PATH, name, fullpath)
        changed = True

    path = op.join(op.expanduser(fullpath), cmd)
    try:
        assert is_exe(path), \
            "***ERROR: Cannot execute binary `{0}`. ".format(path)
    except AssertionError, e:
        if warn == "exit":
            sys.exit("{0!s}Please verify and rerun.".format(e))
        elif warn == "warn":
            logging.warning("{0!s}Some functions may not work.***".format(e))

    if changed:
        configfile = open(cfg, "w")
        config.write(configfile)
        logging.debug("Configuration written to `{0}`.".format(cfg))

    return path


BLPATH = partial(getpath, name="BLAST", url=\
        "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
BDPATH = partial(getpath, name="BEDTOOLS", url=\
        "http://code.google.com/p/bedtools/")
CDPATH = partial(getpath, name="CD-HIT", url=\
        "http://weizhong-lab.ucsd.edu/cd-hit/")
EMBOSSPATH = partial(getpath, name="EMBOSS", url=\
        "http://emboss.sourceforge.net/")
JAVAPATH = partial(getpath, name="Java1.6", url=\
        "http://www.java.com/")
