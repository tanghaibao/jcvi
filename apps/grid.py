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
commitfile = op.join(sge, "COMMIT")
statusfile = op.join(sge, "STATUS")


class CmdSplitter (object):
    """
    The creates parallelized version of cmd, and substituting infile with
    infile_00 and outfile with outfile_00, etc.
    """
    def __init__(self, cmd, N, infile, outfile, outputdir=sge):

        infilepat = re.compile(r"\b%s\b" % infile)
        outfilepat = re.compile(r"\b%s\b" % outfile)

        split_inputs = FileSplitter.get_names(infile, N)
        split_outputs = FileSplitter.get_names(outfile, N)

        # inputs need to be splitted
        fs = FileSplitter(infile, outputdir)
        fs.split(N)

        self.cmds = []

        for sinput, soutput in zip(split_inputs, split_outputs):
            scmd = re.sub(infilepat, sinput, cmd)
            scmd = re.sub(outfilepat, soutput, scmd)
            self.cmds.append(scmd)


class GridProcess (object):

    pat = re.compile(r"Your job (?P<id>[0-9]*) ")

    def __init__(self, cmd, jobid=None):
        self.cmd = cmd
        self.jobid = jobid

    def __str__(self):
        if self.jobid:
            return "\t".join((self.jobid, self.cmd))
        else:
            return self.cmd

    def make_defunct(self):
        self.cmd = "#" + self.cmd

    @property
    def is_defunct(self):
        return self.cmd[0]=='#'

    def start(self):
        # qsub command (the project code is specific to jcvi)
        qsub = "qsub -cwd -P 04048 "

        cmd = qsub + self.cmd
        # run the command and get the job-ID (important)
        p = Popen(cmd, stdout=PIPE, shell=True)
        output = p.communicate()[0]
        
        self.jobid = re.search(self.pat, output).group("id")
        logging.debug("[%s] %s" % (self.jobid, self.cmd))


class Grid (list):

    def __init__(self, cmds=None):

        self.check_sge()

        if cmds:
            for cmd in cmds:
                self.append(GridProcess(cmd))

        else:
            if op.exists(commitfile):
                self.readcommit()
            elif op.exists(statusfile):
                self.readstatus()
            else:
                logging.error("both COMMIT and STATUS not found")
                sys.exit(1)

        assert sum(1 for p in self if not p.is_defunct), "job list need to be non-empty"

    def get_job(self, jobid):
        for p in self:
            if p.jobid==jobid and not p.is_defunct:
                return p

        logging.error("job %s not found in STATUS or is a defunct job" % jobid)
        sys.exit(1)

    def run(self):

        cwd = os.getcwd()

        os.chdir(sge)
        for pi in self:
            pi.start()
        os.chdir(cwd)

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

    def readcommit(self):
        
        if not op.exists(commitfile):
            logging.error("you must run `commit` first")
            sys.exit(1)

        logging.debug("read cmds from %s" % commitfile)
        fp = open(commitfile)
        cmds = [row.strip() for row in fp]
        fp.close()

        for cmd in cmds:
            self.append(GridProcess(cmd))

    def readstatus(self):

        if not op.exists(statusfile):
            logging.error("you mush run `push` first")
            sys.exit(1)

        logging.debug("read cmds from %s" % statusfile)
        fp = open(statusfile)
        cmds = [row.strip().split("\t") for row in fp]
        fp.close()

        for jobid, cmd in cmds:
            self.append(GridProcess(cmd, jobid=jobid))

    def writecommit(self):

        fw = open(commitfile, "w")

        for p in self:
            print >>sys.stderr, p 
            print >>fw, p 

        logging.debug("the above commands are written to %s" % commitfile)
        fw.close()

    def writestatus(self):

        fw = open(statusfile, "w")
        
        for p in self:
            print >>fw, p 
        fw.close()


def main():

    actions = (
        ('commit', 'construct the commands (but do not submit jobs)'),
        ('push', 'run all commands that are commited'),
        ('rerun', 'rerun one command'),
        ('status', 'check status of jobs'),
        ('clean', 'reset the current folder'),
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

    try: # clean up first
        os.remove(commitfile)
        os.remove(statusfile)
    except:
        pass

    cs = CmdSplitter(cmd, N, infile, outfile)

    g = Grid(cs.cmds)
    g.writecommit()


def push(args):
    """
    %prog push

    send the jobs to sun grid engine. the commands are from sge/COMMIT.
    """
    p = OptionParser(push.__doc__)

    g = Grid()
    g.run()
    g.writestatus()

    # remove COMMIT file
    os.remove(commitfile)


def rerun(args):
    """
    %prog rerun jobid
    """
    p = OptionParser(rerun.__doc__)

    p.add_option("-j", dest="jobid", 
            help="rerun job with id (once resubmitted, it will get a new id)")

    opts, args = p.parse_args(args)
    
    if not opts.jobid:
        logging.error("you need to specify jobid (-j)")
        sys.exit(p.print_help())

    g = Grid()
    p = g.get_job(opts.jobid)

    newp = GridProcess(p.cmd)
    newp.start() # start and acquire a job id

    p.make_defunct()
    g.append(newp)
    g.writestatus()


def status(args):
    pass


def clean(args):
    """
    %prog clean

    remove sge folder
    """
    p = OptionParser(clean.__doc__)
    shutil.rmtree(sge)


if __name__ == '__main__':
    main()
