"""
Codes to submit multiple jobs to JCVI grid engine
"""

import os
import os.path as op
import shutil
import sys
import re
import logging
logging.basicConfig(level=logging.DEBUG)

from subprocess import call
from optparse import OptionParser

from jcvi.formats.base import FileSplitter
from jcvi.apps.base import ActionDispatcher

sge = "sge"

class CmdSplitter (object):
    """
    The creates parallelized version of cmd, and substituting infile with
    infile_00 and outfile with outfile_00, etc.
    """
    def __init__(self, cmd, N, infile, outfile):

        infilepat = re.compile(r"\b%s\b" % infile)
        outfilepat = re.compile(r"\b%s\b" % outfile)

        split_inputs = FileSplitter.get_names(infile, N)
        split_outputs = FileSplitter.get_names(outfile, N)

        self.cmds = []

        for sinput, soutput in zip(split_inputs, split_outputs):
            scmd = re.sub(infilepat, sinput, cmd)
            scmd = re.sub(outfilepat, soutput, scmd)
            self.cmds.append(scmd)


class GridProcess (object):

    def __init__(self, cmd):
        self.cmd = cmd

    def __str__(self):
        return self.cmd

    def start(self):
        # qsub command (the project code is specific to jcvi)
        qsub = "qsub -cwd -P 04048 "

        cmd = qsub + self.cmd
        print >>sys.stderr, cmd
        call(cmd, shell=True)


class Grid (list):

    cmds_file = op.join(sge, "cmds")

    def __init__(self, cmds=None):

        self.cmds = cmds or self.readcmds()
        assert self.cmds, "jobs need to be non-empty"

        for cmd in self.cmds:
            self.append(GridProcess(cmd))

    def run(self):
        for pi in self:
            pi.start()

    def check_sge(self, overwrite=False):
        """
        the `sge` folder contains all the jobs
        make sure that this exists; if not create one
        """
        if op.isdir(sge):
            if overwrite:
                shutil.rmtree(sge)
                os.mkdir(sge)
                logging.debug("overwrite %s folder now" % sge)
        else:
            os.mkdir(sge)
            logging.debug("%s folder not found, create one now" % sge)

    def readcmds(self):
        self.check_sge()
        
        logging.debug("read cmds from %s" % self.cmds_file)
        fp = open(self.cmds_file)
        cmds = [row.strip() for row in fp]
        fp.close()

        return cmds

    def writecmds(self):
        self.check_sge()
        self.cmds_file = op.join(sge, "cmds")
        fw = open(self.cmds_file, "w")

        for cmd in self.cmds:
            print >>sys.stderr, cmd
            print >>fw, cmd

        logging.debug("the above commands are written to %s" % \
                self.cmds_file)


def main():

    actions = (
        ('commit', 'construct the commands (but do not submit jobs)'),
        ('push', 'run all commands that are commited'),
        ('rerun', 'rerun one command'),
        ('status', 'check status of jobs')
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def commit(args):

    cmds = CmdSplitter("print infile outfile", 5,
            infile="infile", outfile="outfile").cmds

    g = Grid(cmds)
    g.writecmds()
    pass

def push(args):
    
    g = Grid()
    print g.cmds

def rerun(args):
    pass

def status(args):
    pass



if __name__ == '__main__':
    main()
