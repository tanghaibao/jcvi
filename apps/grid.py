"""
Codes to submit multiple jobs to JCVI grid engine
"""

import os
import os.path as op
import shutil
import sys
import re
import logging

from glob import glob
from itertools import izip_longest
from subprocess import Popen, PIPE
from optparse import OptionParser
from multiprocessing import Process

from jcvi.formats.base import FileSplitter
from jcvi.apps.base import ActionDispatcher, sh, mkdir, debug
debug()


sge = "sge"
PCODE = "04048"  # Project code, JCVI specific
commitfile = op.join(sge, "COMMIT")
statusfile = op.join(sge, "STATUS")


class Jobs (list):
    """
    Runs multiple funcion calls on the SAME computer, using multiprocessing.
    """
    def __init__(self, target, args):

        for x in args:
            self.append(Process(target=target, args=x))

    def run(self):
        for pi in self:
            pi.start()

        for pi in self:
            pi.join()


class CmdReplacer (object):
    """
    Creates parallelized version of cmd, but glob the input files (with the
    same suffix) and produce output files (with another suffix). each command
    runs on one input file and generates corresponding output file
    """
    def __init__(self, cmd, infile, outfile, outputdir=sge):

        split_inputs = glob(infile)
        infile_suffix = op.splitext(infile)[-1]
        outfile_suffix = op.splitext(outfile)[-1]
        split_outputs = [op.basename(x).replace(infile_suffix,
            outfile_suffix) for x in split_inputs]

        self.cmds = []
        for sinput, soutput in zip(split_inputs, split_outputs):
            scmd = cmd.replace(infile, sinput)
            scmd = scmd.replace(outfile, soutput)
            self.cmds.append(scmd)

        assert self.cmds, "no commands generated"

        # cats to different files
        self.outfiles = split_outputs


class CmdSplitter (object):
    """
    This creates parallelized version of cmd, and substituting infile with
    infile_00 and outfile with outfile_00, etc.
    """
    def __init__(self, cmd, infile, outfile, outputdir=sge, N=4):

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

        assert self.cmds, "no commands generated"

        # cats to the same file
        self.outfiles = [outfile] * len(split_outputs)


class GridProcess (object):

    pat = re.compile(r"Your job (?P<id>[0-9]*) ")

    def __init__(self, cmd, jobid="", infile=None, outfile=None, errfile=None):
        self.cmd = cmd
        self.jobid = jobid
        self.infile = infile
        self.outfile = outfile
        self.errfile = errfile

    def __str__(self):
        return "\t".join((x for x in \
                (self.jobid, self.cmd, self.outfile) if x))

    def make_defunct(self):
        if not self.is_defunct:
            self.cmd = "#" + self.cmd

    @property
    def is_defunct(self):
        return self.cmd[0] == '#'

    def start(self, path=sge):

        if self.is_defunct:
            return

        cwd = os.getcwd()
        if path:
            os.chdir(path)

        # qsub command (the project code is specific to jcvi)
        qsub = "qsub -cwd -P {0} ".format(PCODE)

        if self.infile:
            qsub += "-i {0} ".format(self.infile)
        if self.outfile:
            qsub += "-o {0} ".format(self.outfile)
        if self.errfile:
            qsub += "-e {0} ".format(self.errfile)

        cmd = qsub + self.cmd
        # run the command and get the job-ID (important)
        p = Popen(cmd, stdout=PIPE, shell=True)
        output = p.communicate()[0]

        if output.strip() != "":
            self.jobid = re.search(self.pat, output).group("id")
        else:
            self.jobid = "-1"

        msg = "[{0}] {1}".format(self.jobid, self.cmd)
        if self.infile:
            msg += " < {0} ".format(self.infile)
        if self.outfile:
            msg += " > {0} ".format(self.outfile)
        if self.errfile:
            msg += " 2> {0} ".format(self.errfile)

        logging.debug(msg)

        os.chdir(cwd)


class Grid (list):

    def __init__(self, cmds=None, outfiles=[]):

        mkdir(sge)

        if cmds:
            for cmd, outfile in zip(cmds, outfiles):
                self.append(GridProcess(cmd, outfile=outfile))

        else:
            if op.exists(commitfile):
                self.readcommit()
            elif op.exists(statusfile):
                self.readstatus()
            else:
                logging.error("both COMMIT and STATUS not found")
                sys.exit(1)

        assert sum(1 for p in self if not p.is_defunct), \
                "job list need to be non-empty"
        self.cmdgroup = self[0].cmd.split()[0]

    def get_job(self, jobid):
        for p in self:
            if p.jobid == jobid and not p.is_defunct:
                return p

        logging.error("job %s not found in STATUS or is a defunct job" % jobid)
        sys.exit(1)

    def run(self):

        cwd = os.getcwd()

        for pi in self:
            pi.start()

    def readcommit(self):

        if not op.exists(commitfile):
            logging.error("you must run `commit` first")
            sys.exit(1)

        logging.debug("read cmds from %s" % commitfile)
        fp = open(commitfile)
        cmds = [row.strip().split('\t') for row in fp]
        fp.close()

        for cmd, outfile in cmds:
            self.append(GridProcess(cmd, outfile=outfile))

    def readstatus(self):

        if not op.exists(statusfile):
            logging.error("you mush run `push` first")
            sys.exit(1)

        logging.debug("read cmds from %s" % statusfile)
        fp = open(statusfile)
        cmds = [row.strip().split("\t") for row in fp]
        fp.close()

        for jobid, cmd, outfile in cmds:
            self.append(GridProcess(cmd, jobid=jobid, outfile=outfile))

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
        ('run', 'run a normal command on grid'),
        ('commit', 'construct the commands (but do not submit jobs)'),
        ('push', 'run all commands that are commited'),
        ('rerun', 'rerun one command'),
        ('merge', 'merge output files (or stdouts) and stderrs'),
        ('clean', 'reset the current folder'),
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def run(args):
    """
    %prog run command ::: file1 file2

    Parallelize a set of commands on grid. The syntax is modeled after GNU
    parallel <http://www.gnu.org/s/parallel/man.html#options>

    {}   - input line
    {.}  - input line without extension
    {/}  - basename of input line
    {/.} - basename of input line without extension
    {#}  - sequence number of job to run
    :::  - Use arguments from the command line as input source instead of stdin
    (standard input).

    A few examples:
    ls -1 *.fastq | %prog run process {} {.}.pdf  # use stdin
    %prog run process {} {.}.pdf ::: *fastq  # use :::
    %prog run "zcat {} >{.}" ::: *.gz  # quote redirection
    """
    p = OptionParser(run.__doc__)

    if len(args) == 0:
        sys.exit(not p.print_help())

    sep = ":::"
    if sep in args:
        sepidx = args.index(sep)
        filenames = args[sepidx + 1:]
        args = args[:sepidx]
        assert filenames, "No input files"
    else:
        filenames = sys.stdin

    assert args, "Command empty"
    cmd = " ".join(args)

    for i, filename in enumerate(filenames):
        filename = filename.strip()
        noextname = filename.rsplit(".", 1)[0]
        basename = op.basename(filename)
        basenoextname = basename.rsplit(".", 1)[0]
        ncmd = cmd

        if "{}" in ncmd:
            ncmd = ncmd.replace("{}", filename)
        else:
            ncmd += " " + filename

        ncmd = ncmd.replace("{.}", noextname)
        ncmd = ncmd.replace("{/}", basename)
        ncmd = ncmd.replace("{/.}", basenoextname)
        ncmd = ncmd.replace("{#}", str(i))

        outfile = None
        if ">" in ncmd:
            ncmd, outfile = ncmd.split(">", 1)
            ncmd, outfile = ncmd.strip(), outfile.strip()

        p = GridProcess(ncmd, outfile=outfile)
        p.start(path=None)  # current folder


def commit(args):
    """
    %prog commit -i infile -o outfile -n N "command"

    split the command into N pieces, replacing the `infile` string with
    `infile_00` etc., also for `outfile`, write all the commands into
    sge/commit. Most often command needs to be quoted.
    """
    p = OptionParser(commit.__doc__)

    p.add_option("-i", dest="infile", help="infile")
    p.add_option("-o", dest="outfile", help="outfile")
    p.add_option("-n", dest="N", type="int", default=4,
            help="number of hosts to run on (1 < N < 100)")

    opts, args = p.parse_args(args)
    try:
        assert all((opts.infile, opts.outfile)), \
                "need to set all options [-i, -o]"
        assert len(args) == 1
        cmd, N = args[0], opts.N
        infile, outfile = opts.infile, opts.outfile
    except AssertionError as e:
        logging.error(e)
        sys.exit(p.print_help())

    try:  # clean up first
        os.remove(commitfile)
        os.remove(statusfile)
    except:
        pass

    if "*" in opts.infile:
        cs = CmdReplacer(cmd, infile, outfile)
    else:
        cs = CmdSplitter(cmd, infile, outfile, N=N)

    g = Grid(cs.cmds, outfiles=cs.outfiles)
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
    %prog rerun -j jobid
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

    newp = GridProcess(p.cmd, outfile=p.outfile)
    newp.start()  # start and acquire a job id

    p.make_defunct()
    g.append(newp)
    g.writestatus()


def filemerger(split_outputs, outfiles, dir=sge):

    assert len(outfiles) == 1 or len(outfiles) == len(split_outputs), \
            "outfiles must be length 1 or same as split_outputs"

    os.chdir(sge)
    for a, b in izip_longest(split_outputs, outfiles, fillvalue=outfiles[0]):
        sh("cat %s >> %s" % (a, b))


def merge(args):
    """
    %prog merge output|stdout|stderr

    merge all outputfiles into single files, only choices are:
    - output, concatenate all the files that are successful runs in OUTPUT
    - stdout, concatenate all prog.o*
    - stderr, concatenate all prog.e*
    """
    p = OptionParser(merge.__doc__)

    opts, args = p.parse_args(args)
    try:
        filetype = args[0]
        assert filetype in ("output", "stdout", "stderr")
    except Exception, e:
        logging.error(e)
        sys.exit(p.print_help())

    assert op.exists(statusfile), "STATUS file not found, cannot merge"

    g = Grid()
    cmdgroup = g.cmdgroup  # the actual command that was run
    stdoutlist = []
    stderrlist = []
    outfiles = []
    for p in g:
        if p.is_defunct:
            continue
        stdoutlist.append("%s.o%s" % (cmdgroup, p.jobid))
        stderrlist.append("%s.e%s" % (cmdgroup, p.jobid))
        outfiles.append(p.outfile)

    if filetype == "output":
        filemerger(stdoutlist, outfiles)
    elif filetype == "stdout":
        filemerger(stdoutlist, ["%s.stdout" % cmdgroup])
    elif filetype == "stderr":
        filemerger(stderrlist, ["%s.stderr" % cmdgroup])


def clean(args):
    """
    %prog clean

    remove sge folder
    """
    p = OptionParser(clean.__doc__)

    if op.exists(sge):
        shutil.rmtree(sge)


if __name__ == '__main__':
    main()
