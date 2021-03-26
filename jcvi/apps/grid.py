"""
Codes to submit multiple jobs to JCVI grid engine
"""
from __future__ import print_function

import os.path as op
import sys
import re
import logging
import platform

from multiprocessing import (
    Pool,
    Process,
    Value,
    cpu_count,
    get_context,
    set_start_method,
)
from multiprocessing.queues import Queue

from jcvi.formats.base import write_file, must_open
from jcvi.apps.base import (
    OptionParser,
    ActionDispatcher,
    popen,
    backup,
    mkdir,
    sh,
    listify,
)


class SharedCounter(object):
    """A synchronized shared counter.

    The locking done by multiprocessing.Value ensures that only a single
    process or thread may read or write the in-memory ctypes object. However,
    in order to do n += 1, Python performs a read followed by a write, so a
    second process may read the old value before the new one is written by the
    first process. The solution is to use a multiprocessing.Lock to guarantee
    the atomicity of the modifications to Value.

    This class comes almost entirely from Eli Bendersky's blog:
    http://eli.thegreenplace.net/2012/01/04/shared-counter-with-pythons-multiprocessing/
    """

    def __init__(self, n=0):
        self.count = Value("i", n)

    def increment(self, n=1):
        """ Increment the counter by n (default = 1) """
        with self.count.get_lock():
            self.count.value += n

    @property
    def value(self):
        """ Return the value of the counter """
        return self.count.value


class Queue(Queue):
    """A portable implementation of multiprocessing.Queue.

    Because of multithreading / multiprocessing semantics, Queue.qsize() may
    raise the NotImplementedError exception on Unix platforms like Mac OS X
    where sem_getvalue() is not implemented. This subclass addresses this
    problem by using a synchronized shared counter (initialized to zero) and
    increasing / decreasing its value every time the put() and get() methods
    are called, respectively. This not only prevents NotImplementedError from
    being raised, but also allows us to implement a reliable version of both
    qsize() and empty().
    """

    def __init__(self, *args, **kwargs):
        super(Queue, self).__init__(*args, **kwargs, ctx=get_context())
        self.size = SharedCounter(0)

    def put(self, *args, **kwargs):
        self.size.increment(1)
        super(Queue, self).put(*args, **kwargs)

    def get(self, *args, **kwargs):
        self.size.increment(-1)
        return super(Queue, self).get(*args, **kwargs)

    def qsize(self):
        """ Reliable implementation of multiprocessing.Queue.qsize() """
        return self.size.value

    def empty(self):
        """ Reliable implementation of multiprocessing.Queue.empty() """
        return not self.qsize()


class Parallel(object):
    """
    Run a number of commands in parallel.
    """

    def __init__(self, cmds, cpus=cpu_count()):
        self.cmds = cmds
        self.cpus = min(len(cmds), cpus)

    def run(self):
        p = Pool(processes=self.cpus)
        p.map(sh, self.cmds)


class Dependency(object):
    """
    Used by MakeManager.
    """

    def __init__(self, source, target, cmds, id, remove=False):
        self.id = id
        self.source = listify(source)
        self.target = listify(target)
        self.cmds = listify(cmds)
        if remove:
            rm_cmd = "rm -f {0}".format(" ".join(self.target))
            self.cmds = [rm_cmd] + self.cmds

    def __str__(self):
        source = " ".join(self.source)
        target = " ".join(self.target)
        # When there are multiple targets, use .INTERMEDIATE
        # <http://stackoverflow.com/questions/2973445/gnu-makefile-rule-generating-a-few-targets-from-a-single-source-file>
        if len(self.target) > 1:
            intermediate = "{0}.intermediate".format(self.id)
            s = "{0} : {1}\n".format(target, intermediate)
            s += ".INTERMEDIATE: {0}\n".format(intermediate)
            s += "{0} : {1}\n".format(intermediate, source)
        else:
            s = "{0} : {1}\n".format(target, source)

        for c in self.cmds:
            c = c.replace("$", "$$")  # Command escaping
            s += "\t" + c + "\n"
        return s


class MakeManager(list):
    """
    Write and execute makefile.
    """

    def __init__(self, filename="makefile"):
        self.makefile = filename
        self.targets = set()
        self.ndeps = 0

    def add(self, source, target, cmds, remove=False):
        self.ndeps += 1
        d = Dependency(source, target, cmds, self.ndeps, remove=remove)
        self.append(d)
        self.targets |= set(listify(target))

    def write(self):
        assert self.targets, "No targets specified"
        filename = self.makefile
        if op.exists(filename):
            backup(filename)
        fw = open(filename, "w")
        print("all : {0}\n".format(" ".join(sorted(self.targets))), file=fw)
        for d in self:
            print(d, file=fw)
        print("clean :\n\trm -rf {0}\n".format(" ".join(self.targets)), file=fw)
        fw.close()
        logging.debug("Makefile written to `{0}`.".format(self.makefile))

    def run(self, cpus=1):
        if not op.exists(self.makefile):
            self.write()
        cmd = "make -j {0} -f {1}".format(cpus, self.makefile)
        sh(cmd)

    def clean(self):
        cmd = "make clean -f {}".format(self.makefile)
        sh(cmd)


class Jobs(list):
    """
    Runs multiple funcion calls on the SAME computer, using multiprocessing.
    """

    def __init__(self, target, args):

        for x in args:
            x = listify(x)
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


class Poison:
    pass


class WriteJobs(object):
    """
    Runs multiple function calls, but write to the same file.

    Producer-consumer model.
    """

    def __init__(self, target, args, filename, cpus=cpu_count()):
        # macOS starts process with fork by default: https://zhuanlan.zhihu.com/p/144771768
        if platform.system() == "Darwin":
            set_start_method("fork")

        workerq = Queue()
        writerq = Queue()

        for a in args:
            workerq.put(a)

        cpus = min(cpus, len(args))
        for i in range(cpus):
            workerq.put(Poison())

        self.worker = Jobs(work, args=[(workerq, writerq, target)] * cpus)
        self.writer = Process(target=write, args=(workerq, writerq, filename, cpus))

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
    from rich.progress import Progress

    fw = must_open(filename, "w")
    isize = queue_in.qsize()
    logging.debug("A total of {0} items to compute.".format(isize))
    isize = isize or 1
    poisons = 0
    with Progress() as progress:
        task = progress.add_task("[green]Processing ...", total=isize)
        while True:
            res = queue_out.get()
            qsize = queue_in.qsize()
            progress.update(task, completed=isize - qsize)
            if isinstance(res, Poison):
                poisons += 1
                if poisons == cpus:  # wait all workers finish
                    break
            elif res:
                print(res, file=fw)
                fw.flush()
    fw.close()


class GridOpts(dict):
    def __init__(self, opts):
        export = (
            "pcode",
            "queue",
            "threaded",
            "concurrency",
            "outdir",
            "name",
            "hold_jid",
        )
        for e in export:
            if e in opts.__dict__:
                self[e] = getattr(opts, e)


class GridProcess(object):

    pat1 = re.compile(r"Your job (?P<id>[0-9]*) ")
    pat2 = re.compile(r"Your job-array (?P<id>\S*) ")

    def __init__(
        self,
        cmd,
        jobid="",
        pcode="99999",
        queue="default",
        threaded=None,
        infile=None,
        outfile=None,
        errfile=None,
        arr=None,
        concurrency=None,
        outdir=".",
        name=None,
        hold_jid=None,
        extra_opts=None,
        grid_opts=None,
    ):

        self.cmd = cmd
        self.jobid = jobid
        self.queue = queue
        self.threaded = threaded
        self.infile = infile
        self.outfile = outfile or ""
        self.errfile = errfile or ""
        self.arr = arr
        self.concurrency = concurrency
        self.outdir = outdir
        self.name = name
        self.pcode = pcode
        self.hold_jid = hold_jid
        self.pat = self.pat2 if arr else self.pat1
        self.extra = extra_opts if extra_opts else None
        if grid_opts:
            self.__dict__.update(GridOpts(grid_opts))

    def __str__(self):
        return "\t".join((x for x in (self.jobid, self.cmd, self.outfile) if x))

    def build(self):
        # Shell commands
        if "|" in self.cmd or "&&" in self.cmd or "||" in self.cmd:
            quote = '"' if "'" in self.cmd else "'"
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
        if self.hold_jid:
            param = "-hold_jid_ad" if self.arr else "-hold_jid"
            qsub += " {0} {1}".format(param, self.hold_jid)
        if self.extra:
            qsub += " {0}".format(self.extra)

        # I/O
        infile = self.infile
        outfile = self.outfile
        errfile = self.errfile
        outdir = self.outdir
        mkdir(outdir)
        redirect_same = outfile and (outfile == errfile)

        if infile:
            qsub += " -i {0}".format(infile)
        if outfile:
            self.outfile = op.join(outdir, outfile)
            qsub += " -o {0}".format(self.outfile)
        if errfile:
            if redirect_same:
                qsub += " -j y"
            else:
                self.errfile = op.join(outdir, errfile)
                qsub += " -e {0}".format(self.errfile)

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
            msg += " > {0} ".format(self.outfile)
        if self.errfile:
            backup(self.errfile)
            msg += " 2> {0} ".format(self.errfile)

        logging.debug(msg)


class Grid(list):
    def __init__(self, cmds, outfiles=[]):

        assert cmds, "Commands empty!"
        if not outfiles:
            outfiles = [None] * len(cmds)

        for cmd, outfile in zip(cmds, outfiles):
            self.append(GridProcess(cmd, outfile=outfile))

    def run(self):
        for pi in self:
            pi.start()


PBS_STANZA = """
#PBS -q standard
#PBS -J 1-{0}
#PBS -l select=1:ncpus={1}:mem=23gb
#PBS -l pvmem=23gb
#PBS -l walltime=100:00:00
#PBS -W group_list=genomeanalytics
"""

arraysh = """
CMD=`awk "NR==$SGE_TASK_ID" {0}`
$CMD"""

arraysh_ua = (
    PBS_STANZA
    + """
cd $PBS_O_WORKDIR
CMD=`awk "NR==$PBS_ARRAY_INDEX" {2}`
$CMD"""
)


def get_grid_engine():
    cmd = "qsub --version"
    ret = popen(cmd, debug=False).read()
    return "PBS" if "PBS" in ret else "SGE"


def main():

    actions = (
        ("run", "run a normal command on grid"),
        ("array", "run an array job"),
        ("kill", "wrapper around the `qdel` command"),
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
    p.set_params(prog="grid")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (cmds,) = args
    fp = open(cmds)
    N = sum(1 for x in fp)
    fp.close()

    pf = cmds.rsplit(".", 1)[0]
    runfile = pf + ".sh"
    assert runfile != cmds, "Commands list file should not have a `.sh` extension"

    engine = get_grid_engine()
    threaded = opts.threaded or 1
    contents = (
        arraysh.format(cmds)
        if engine == "SGE"
        else arraysh_ua.format(N, threaded, cmds)
    )
    write_file(runfile, contents)

    if engine == "PBS":
        return

    outfile = "{0}.{1}.out".format(pf, "\$TASK_ID")
    errfile = "{0}.{1}.err".format(pf, "\$TASK_ID")
    p = GridProcess(
        "sh {0}".format(runfile),
        outfile=outfile,
        errfile=errfile,
        arr=N,
        extra_opts=opts.extra,
        grid_opts=opts,
    )
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
    p.set_params(prog="grid")
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    sep = ":::"
    if sep in args:
        sepidx = args.index(sep)
        filenames = args[sepidx + 1 :]
        args = args[:sepidx]
        if not filenames:
            filenames = [""]
    else:
        filenames = sys.stdin if not sys.stdin.isatty() else [""]

    cmd = " ".join(args)

    cmds = [] if filenames else [(cmd, None)]
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
        cmds.append((ncmd, outfile))

    for ncmd, outfile in cmds:
        p = GridProcess(ncmd, outfile=outfile, extra_opts=opts.extra, grid_opts=opts)
        p.start()


def guess_method(tag):
    from jcvi.formats.base import is_number

    jobids = tag.split(",")
    for jobid in jobids:
        if not is_number(jobid):
            return "pattern"
    return "jobid"


def kill(args):
    """
    %prog kill [options] JOBNAMEPAT/JOBIDs

    Kill jobs based on JOBNAME pattern matching (case-sensitive)
    or list of JOBIDs (comma separated)

    Examples:
    %prog kill "pyth*"                 # Use regex
    %prog kill 160253,160245,160252    # Use list of job ids
    %prog kill all                     # Everything
    """
    import shlex
    from jcvi.apps.base import sh, getusername
    from subprocess import check_output, CalledProcessError
    import xml.etree.ElementTree as ET

    valid_methods = ("pattern", "jobid")
    p = OptionParser(kill.__doc__)
    p.add_option(
        "--method",
        choices=valid_methods,
        help="Identify jobs based on [default: guess]",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    username = getusername()
    (tag,) = args
    tag = tag.strip()

    if tag == "all":
        sh("qdel -u {0}".format(username))
        return

    valid_jobids = set()
    method = opts.method or guess_method(tag)
    if method == "jobid":
        jobids = tag.split(",")
        valid_jobids |= set(jobids)
    elif method == "pattern":
        qsxmlcmd = 'qstat -u "{}" -j "{}" -nenv -njd -xml'.format(username, tag)
        try:
            qsxml = check_output(shlex.split(qsxmlcmd)).strip()
        except CalledProcessError as e:
            qsxml = None
            logging.debug('No jobs matching the pattern "{0}"'.format(tag))

        if qsxml is not None:
            for job in ET.fromstring(qsxml).findall("djob_info"):
                for elem in job.findall("element"):
                    jobid = elem.find("JB_job_number").text
                    valid_jobids.add(jobid)

    if valid_jobids:
        sh("qdel {0}".format(",".join(valid_jobids)))


if __name__ == "__main__":
    main()
