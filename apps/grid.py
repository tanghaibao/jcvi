"""
Codes to submit multiple jobs to JCVI grid engine
"""

import os
import os.path as op
import sys
import time
import re
import logging

from subprocess import Popen, PIPE
from jcvi.apps.base import OptionParser
from multiprocessing import Process, Queue, cpu_count

from jcvi.formats.base import write_file, must_open
from jcvi.apps.base import ActionDispatcher, sh, popen, backup, mkdir, debug
debug()



class Jobs (list):
    """
    Runs multiple funcion calls on the SAME computer, using multiprocessing.
    """
    def __init__(self, target, args):

        for x in args:
            self.append(Process(target=target, args=x))

    def start(self):
        for pi in self:
            pi.start()

    def join(self):
        for pi in self:
            pi.join()

    def run(self):
        self.start()
        self.join()


class Poison (int):
    def __init__(self):
        pass


class WriteJobs (object):
    """
    Runs multiple function calls, but write to the same file.

    Producer-consumer model.
    """
    def __init__(self, target, args, filename, cpus=cpu_count()):
        workerq = Queue()
        writerq = Queue()

        for a in args:
            workerq.put(a)
        logging.debug("A total of {0} items to compute.".format(workerq.qsize()))

        for i in xrange(cpus):
            workerq.put(Poison())

        self.worker = Jobs(work, args=[(workerq, writerq, target)] * cpus)
        self.writer = Process(target=write, args=(workerq, writerq, \
                                                  filename, cpus))

    def run(self):
        self.worker.start()
        self.writer.start()
        self.worker.join()
        self.writer.join()


def work(queue_in, queue_out, target):
    while True:
        a = queue_in.get()
        if isinstance(a, Poison):
            break
        res = target(a)
        queue_out.put(res)
    queue_out.put(Poison())


def write(queue_in, queue_out, filename, cpus):
    fw = must_open(filename, "w")
    poisons = 0
    while True:
        res = queue_out.get()
        logging.debug("Writer: {0} items left to compute".\
                        format(queue_in.qsize()))
        if isinstance(res, Poison):
            poisons += 1
            if poisons == cpus:  # wait all workers finish
                break
        if res:
            print >> fw, res
            fw.flush()
    fw.close()


class GridProcess (object):

    pat1 = re.compile(r"Your job (?P<id>[0-9]*) ")
    pat2 = re.compile(r"Your job-array (?P<id>\S*) ")

    def __init__(self, cmd, jobid="", pcode="04048", queue="default", threaded=None,
                       infile=None, outfile=None, errfile=None, arr=None,
                       concurrency=None, outdir=".", name=None):

        self.cmd = cmd
        self.jobid = jobid
        self.queue = queue
        self.threaded = threaded
        self.infile = infile
        self.outfile = outfile
        self.errfile = errfile
        self.arr = arr
        self.concurrency = concurrency
        self.outdir = outdir
        self.name = name
        self.pcode = pcode
        self.pat = self.pat2 if arr else self.pat1

    def __str__(self):
        return "\t".join((x for x in \
                (self.jobid, self.cmd, self.outfile) if x))

    def build(self):
        # Shell commands
        if "|" in self.cmd or "&&" in self.cmd or "||" in self.cmd:
            quote = "\"" if "'" in self.cmd else "'"
            self.cmd = "sh -c {1}{0}{1}".format(self.cmd, quote)

        # qsub command (the project code is specific to jcvi)
        qsub = "qsub -P {0} -cwd".format(self.pcode)
        if self.queue != "default":
            qsub += " -l {0}".format(self.queue)
        if self.threaded:
            qsub += " -pe threaded {0}".format(self.threaded)
        if self.arr:
            assert 1 <= self.arr < 100000
            qsub += " -t 1-{0}".format(self.arr)
        if self.concurrency:
            qsub += " -tc {0}".format(self.concurrency)
        if self.name:
            qsub += ' -N "{0}"'.format(self.name)

        # I/O
        infile = self.infile
        outfile = self.outfile
        errfile = self.errfile
        outdir = self.outdir
        mkdir(outdir)
        redirect_same = outfile and (outfile == errfile)

        if infile:
            qsub += " -i {0}".format(infile)
        if redirect_same:
            qsub += " -j y"
        if outfile:
            outfile = op.join(outdir, outfile)
            qsub += " -o {0}".format(outfile)
        if errfile and not redirect_same:
            errfile = op.join(outdir, errfile)
            qsub += " -e {0}".format(errfile)

        cmd = " ".join((qsub, self.cmd))
        return cmd

    def start(self):
        cmd = self.build()
        # run the command and get the job-ID (important)
        output = popen(cmd, debug=False).read()

        if output.strip() != "":
            self.jobid = re.search(self.pat, output).group("id")
        else:
            self.jobid = "-1"

        msg = "[{0}] {1}".format(self.jobid, self.cmd)
        if self.infile:
            msg += " < {0} ".format(self.infile)
        if self.outfile:
            backup(self.outfile)
            msg += " > {0}/{1} ".format(self.outdir, self.outfile)
        if self.errfile:
            backup(self.errfile)
            msg += " 2> {0}/{1} ".format(self.outdir, self.errfile)

        logging.debug(msg)


class Grid (list):

    def __init__(self, cmds, outfiles=[]):

        assert cmds, "Commands empty!"
        if not outfiles:
            outfiles = [None] * len(cmds)

        for cmd, outfile in zip(cmds, outfiles):
            self.append(GridProcess(cmd, outfile=outfile))

    def run(self):

        cwd = os.getcwd()

        for pi in self:
            pi.start()


arraysh = """
CMD=`awk "NR==$SGE_TASK_ID" {0}`
$CMD"""


def main():

    actions = (
        ('run', 'run a normal command on grid'),
        ('array', 'run an array job'),
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def array(args):
    """
    %prog array commands.list

    Parallelize a set of commands on grid using array jobs.
    """
    p = OptionParser(array.__doc__)
    p.set_grid_opts(array=True)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    cmds, = args
    fp = open(cmds)
    ncmds = sum(1 for x in fp)
    fp.close()

    pf = cmds.rsplit(".",  1)[0]
    runfile = pf + ".sh"
    assert runfile != cmds, \
            "Commands list file should not have a `.sh` extension"

    contents = arraysh.format(cmds)
    write_file(runfile, contents, meta="run script")

    outfile = "{0}.{1}.out".format(pf, "\$TASK_ID")
    p = GridProcess("sh {0}".format(runfile), outfile=outfile, errfile=outfile,
                    pcode=opts.pcode, queue=opts.queue, threaded=opts.threaded,
                    arr=ncmds, concurrency=opts.concurrency,
                    outdir=opts.outdir, name=opts.name)
    p.start()


def run(args):
    """
    %prog run command ::: file1 file2

    Parallelize a set of commands on grid. The syntax is modeled after GNU
    parallel <http://www.gnu.org/s/parallel/man.html#options>

    {}   - input line
    {.}  - input line without extension
    {_}  - input line first part
    {/}  - basename of input line
    {/.} - basename of input line without extension
    {/_} - basename of input line first part
    {#}  - sequence number of job to run
    :::  - Use arguments from the command line as input source instead of stdin
    (standard input).

    If file name is `t/example.tar.gz`, then,
    {} is "t/example.tar.gz", {.} is "t/example.tar", {_} is "t/example"
    {/} is "example.tar.gz", {/.} is "example.tar", {/_} is "example"

    A few examples:
    ls -1 *.fastq | %prog run process {} {.}.pdf  # use stdin
    %prog run process {} {.}.pdf ::: *fastq  # use :::
    %prog run "zcat {} > {.}" ::: *.gz  # quote redirection
    %prog run < commands.list  # run a list of commands
    """
    p = OptionParser(run.__doc__)
    p.set_grid_opts()
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    sep = ":::"
    if sep in args:
        sepidx = args.index(sep)
        filenames = args[sepidx + 1:]
        args = args[:sepidx]
        if not filenames:
            filenames = [""]
    else:
        filenames = sys.stdin

    cmd = " ".join(args)

    for i, filename in enumerate(filenames):
        filename = filename.strip()
        noextname = filename.rsplit(".", 1)[0]
        prefix, basename = op.split(filename)
        basenoextname = basename.rsplit(".", 1)[0]
        basefirstname = basename.split(".")[0]
        firstname = op.join(prefix, basefirstname)
        ncmd = cmd

        if "{" in ncmd:
            ncmd = ncmd.replace("{}", filename)
        else:
            ncmd += " " + filename

        ncmd = ncmd.replace("{.}", noextname)
        ncmd = ncmd.replace("{_}", firstname)
        ncmd = ncmd.replace("{/}", basename)
        ncmd = ncmd.replace("{/.}", basenoextname)
        ncmd = ncmd.replace("{/_}", basefirstname)
        ncmd = ncmd.replace("{#}", str(i))

        outfile = None
        if ">" in ncmd:
            ncmd, outfile = ncmd.split(">", 1)
            ncmd, outfile = ncmd.strip(), outfile.strip()

        ncmd = ncmd.strip()
        p = GridProcess(ncmd, outfile=outfile, pcode=opts.pcode,
                        queue=opts.queue, threaded=opts.threaded,
                        outdir=opts.outdir, name=opts.name)
        p.start()


if __name__ == '__main__':
    main()
