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

from subprocess import Popen, PIPE 
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

    pat = re.compile(r"Your job (?P<id>[0-9]*) ")

    def __init__(self, cmd):
        self.cmd = cmd

    def __str__(self):
        return self.cmd

    def start(self):
        # qsub command (the project code is specific to jcvi)
        qsub = "qsub -cwd -P 04048 "

        cmd = qsub + self.cmd
        # run the command and get the job-ID (important)
        p = Popen(cmd, stdout=PIPE)
        output = p.communicate()[0]
        
        self.jobid = re.search(self.pat, output).group("id")
        logging.debug("[%s] %s" % (self.jobid, self.cmd))


class Grid (list):

    commitfile = op.join(sge, "COMMIT")
    statusfile = op.join(sge, "STATUS")

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
        fp = open(self.commitfile)
        cmds = [row.strip() for row in fp]
        fp.close()

        return cmds

    def writecmds(self):
        self.check_sge()
        fw = open(self.commitfile, "w")

        for cmd in self.cmds:
            print >>sys.stderr, cmd
            print >>fw, cmd

        logging.debug("the above commands are written to %s" % \
                self.commitfile)

    def writejobids(self):
        self.check_sge()
        fw = open(self.statusfile, "w")
        
        for ps in self:
            print >>fw, "%s %s" % (p.jobid, p.cmd)


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
    """
    %prog commit -i infile -o outfile -c "command" -n N

    split the command into N pieces, replacing the `infile` string with
    `infile_00` etc., also for `outfile`, write all the commands into
    sge/commit
    """
    p = OptionParser(commit.__doc__)
    
    p.add_option("-i", dest="infile", help="infile")
    p.add_option("-o", dest="outfile", help="outfile")
    p.add_option("-c", dest="cmd", 
            help="cmd needs to parallelize (please quote the string)")
    p.add_option("-n", dest="N", type="int", 
            help="number of hosts to run on (1 < N < 100)")

    opts, args = p.parse_args(args)
    try:
        assert all((opts.cmd, opts.N, opts.infile, opts.outfile)), \
                "need to set all options [-i, -o, -c, -n]"
        cmd, N = opts.cmd, opts.N 
        infile, outfile = opts.infile, opts.outfile
    except Exception, e:
        logging.error(str(e))
        sys.exit(p.print_help())

    cs = CmdSplitter(cmd, N, infile, outfile)

    g = Grid(cs.cmds)
    g.writecmds()


def push(args):
    """
    %prog push

    send the jobs to sun grid engine. the commands are from sge/COMMIT.
    """
    g = Grid()
    g.run()


def rerun(args):
    """
    %prog rerun jobid
    """
    p = OptionParser(rerun.__doc__)

    p.add_option("-j", dest="job", type="int",
            help="rerun job with id (once resubmitted, it will get a new id)")

    g = Grid()
    g[i].start()


def status(args):
    pass


if __name__ == '__main__':
    main()
