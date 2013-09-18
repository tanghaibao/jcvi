"""
basic support for running library as script
"""

import os
import os.path as op
import sys
import shutil
import logging

from subprocess import call
from optparse import OptionParser as OptionP, OptionGroup

os.environ["LC_ALL"] = "C"


class ActionDispatcher (object):
    """
    This class will be invoked
    a) when either a directory is run via __main__, listing all SCRIPTs
    b) when a script is run directly, listing all ACTIONs

    This is controlled through the meta variable, which is automatically
    determined in get_meta().
    """
    def __init__(self, actions):

        self.actions = actions
        if not actions:
            actions = [(None, None)]
        self.valid_actions, self.action_helps = zip(*actions)

    def get_meta(self):
        args = splitall(sys.argv[0])[-3:]
        args[-1] = args[-1].replace(".py", "")
        meta = "SCRIPT" if args[-1] == "__main__" else "ACTION"
        return meta, args

    def print_help(self):
        meta, args = self.get_meta()
        if meta == "SCRIPT":
            args[-1] = meta
        else:
            args[-1] += " " + meta

        help = "Usage:\n    python -m {0}\n\n".format(".".join(args))
        help += "Available {0}s:\n".format(meta)
        for action, action_help in self.actions:
            help += "    `%s`: %s\n" % (action, action_help)

        print >> sys.stderr, help
        sys.exit(1)

    def dispatch(self, globals):
        meta = "ACTION"  # function is only invoked for listing ACTIONs
        if len(sys.argv) == 1:
            self.print_help()

        action = sys.argv[1]

        if not action in self.valid_actions:
            print >> sys.stderr, "[error] {0} not a valid {1}\n".format(action, meta)
            self.print_help()

        globals[action](sys.argv[2:])


class OptionParser (OptionP):

    def __init__(self, doc):

        OptionP.__init__(self, doc)

    def set_grid(self):
        """
        Add --grid options for command line programs
        """
        self.add_option("--grid", dest="grid",
                default=False, action="store_true",
                help="Run on the grid [default: %default]")

    def set_grid_opts(self, array=False):
        queue_choices = ("default", "fast", "medium", "himem")
        valid_pcodes = popen("qconf -sprjl", debug=False).read().strip().split("\n")
        self.add_option("-P", dest="pcode", default="04048", choices=valid_pcodes,
                     help="Specify accounting project code [default: %default]")
        self.add_option("-l", dest="queue", default="default", choices=queue_choices,
                     help="Name of the queue, one of {0} [default: %default]". \
                          format("|".join(queue_choices)))
        self.add_option("-t", dest="threaded", type="int",
                     help="Append '-pe threaded N' [default: %default]")
        if array:
            self.add_option("-c", dest="concurrency", type="int",
                     help="Append task concurrency limit '-tc N' [default: %default]")
        self.add_option("-d", dest="outdir", default=".",
                help="Specify directory to store grid output/error files" +
                      " [default: %default]")
        self.add_option("-N", dest="name", default=None,
                     help="Specify descriptive name for the job [default: %default]")

    def set_params(self):
        """
        Add --params options for given command line programs
        """
        self.add_option("--params", dest="extra", default="",
                help="Extra parameters to run [default: %default]")

    def set_outfile(self, outfile="stdout"):
        """
        Add --outfile options to print out to filename.
        """
        self.add_option("-o", "--outfile", default=outfile,
                help="Outfile name [default: %default]")

    def set_tmpdir(self, tmpdir="/tmp"):
        """
        Add --temporary_directory option to specify unix `sort` tmpdir
        """
        self.add_option("-T", "--tmpdir", default=tmpdir,
                help="Use temp directory instead of $TMP [default: %default]")

    def set_cpus(self, cpus=0):
        """
        Add --cpus options to specify how many threads to use.
        """
        from multiprocessing import cpu_count

        max_cpus = cpu_count()
        if not 0 < cpus < max_cpus:
            cpus = max_cpus
        self.add_option("--cpus", default=cpus, type="int",
                     help="Number of CPUs to use, 0=unlimited [default: %default]")

    def set_stripnames(self):
        self.add_option("--no_strip_names", dest="strip_names",
                action="store_false", default=True,
                help="do not strip alternative splicing "
                "(e.g. At5g06540.1 -> At5g06540)")

    def set_beds(self):
        self.add_option("--qbed", help="Path to qbed")
        self.add_option("--sbed", help="Path to sbed")

    def set_sam_options(self):
        self.add_option("--bam", default=False, action="store_true",
                     help="write to bam file [default: %default]")
        self.add_option("--uniq", default=False, action="store_true",
                     help="Keep only uniquely mapped [default: %default]")
        self.add_option("--unmapped", default=False, action="store_true",
                     help="Keep unmapped reads [default: %default]")
        self.set_cpus()
        self.set_params()
        self.set_grid()

    def set_min_length(self, default=10000):
        self.add_option("--min_length", default=default, type="int",
                     help="Minimum length to consider [default: %default]")

    def set_mingap(self, default=100):
        self.add_option("--mingap", default=default, type="int",
                     help="Minimum size of gaps [default: %default]")

    def set_align(self, pctid="off", hitlen="off", pctcov="off", evalue="off"):
        if pctid != "off":
            self.add_option("--pctid", default=pctid, type="int",
                     help="Sequence percent identity [default: %default]")
        if hitlen != "off":
            self.add_option("--hitlen", default=hitlen, type="int",
                     help="Minimum overlap length [default: %default]")
        if pctcov != "off":
            self.add_option("--pctcov", default=pctcov, type="int",
                     help="Percentage coverage cutoff [default: %default]")
        if evalue != "off":
            self.add_option("--evalue", default=evalue, type="float",
                     help="E-value cutoff [default: %default]")


def ConfigSectionMap(Config, section):
    """
    Read a specific section from a ConfigParser() object and return
    a dict() of all key-value pairs in that section
    """
    debug()

    cfg = {}
    options = Config.options(section)
    for option in options:
        try:
            cfg[option] = Config.get(section, option)
            if cfg[option] == -1:
                logging.debug("skip: {0}".format(option))
        except:
            logging.debug("exception on {0}!".format(option))
            cfg[option] = None
    return cfg


def splitall(path):
    allparts = []
    while True:
        path, p1 = op.split(path)
        if not p1:
            break
        allparts.append(p1)
    allparts = allparts[::-1]
    return allparts


def get_module_docstring(filepath):
    "Get module-level docstring of Python module at filepath, e.g. 'path/to/file.py'."
    co = compile(open(filepath).read(), filepath, 'exec')
    if co.co_consts and isinstance(co.co_consts[0], basestring):
        docstring = co.co_consts[0]
    else:
        docstring = None
    return docstring


def dmain(cwd):
    pyscripts = glob(op.join(cwd, "*.py"))
    actions = []
    for ps in sorted(pyscripts):
        action = op.basename(ps).replace(".py", "")
        if action[0] == "_":  # hidden namespace
            continue
        pd = get_module_docstring(ps)
        action_help = [x.rstrip(":.,\n") for x in pd.splitlines(True) \
                if len(x.strip()) > 10][0] if pd else "no docstring found"
        actions.append((action, action_help))

    a = ActionDispatcher(actions)
    a.print_help()


def backup(filename):
    if op.exists(filename):
        bakname = filename + ".bak"
        logging.debug("Backup `{0}` to `{1}`".format(filename, bakname))
        sh("mv {0} {1}".format(filename, bakname))
        return bakname


def gethostname():
    from socket import gethostname
    return gethostname()


def getusername():
    from getpass import getuser
    return getuser()


def getdomainname():
    from socket import getfqdn
    return ".".join(str(x) for x in getfqdn().split(".")[1:])


def sh(cmd, grid=False, infile=None, outfile=None, errfile=None,
        append=False, background=False, threaded=None, log=True,
        shell="/bin/bash"):
    """
    simple wrapper for system calls
    """
    if grid:
        from jcvi.apps.grid import GridProcess
        pr = GridProcess(cmd, infile=infile, outfile=outfile, errfile=errfile,
                         threaded=threaded)
        pr.start()
        return 0  # A fake retcode
    else:
        if infile:
            cat = "cat"
            if infile.endswith(".gz"):
                cat = "zcat"
            cmd = "{0} {1} | ".format(cat, infile) + cmd
        if outfile and outfile != "stdout":
            if outfile.endswith(".gz"):
                cmd += " | gzip"
            tag = ">"
            if append:
                tag = ">>"
            cmd += " {0} {1} ".format(tag, outfile)
        if errfile:
            cmd += " 2> {0} ".format(errfile)
        if background:
            cmd += " & "

        if log:
            logging.debug(cmd)
        return call(cmd, shell=True, executable=shell)


def popen(cmd, debug=True, shell="/bin/bash"):
    """
    Capture the cmd stdout output to a file handle.
    """
    from subprocess import Popen, PIPE
    if debug:
        logging.debug(cmd)
    proc = Popen(cmd, bufsize=1, stdout=PIPE, shell=True, \
                 executable=shell)
    return proc.stdout


def is_exe(fpath):
    return op.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """
    Emulates the unix which command.

    >>> which("cat")
    "/bin/cat"
    >>> which("nosuchprogram")
    """
    fpath, fname = op.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = op.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def glob(pathname):
    """
    Wraps around glob.glob(), but return a sorted list.
    """
    import glob as gl
    return sorted(gl.glob(pathname))


def mkdir(dirname, overwrite=False):
    """
    Wraps around os.mkdir(), but checks for existence first.
    """
    if op.isdir(dirname):
        if overwrite:
            shutil.rmtree(dirname)
            os.mkdir(dirname)
            logging.debug("Overwrite folder `{0}`.".format(dirname))
        else:
            return False  # Nothing is changed
    else:
        try:
            os.mkdir(dirname)
        except:
            os.makedirs(dirname)
        logging.debug("`{0}` not found. Creating new.".format(dirname))

    return True


def is_newer_file(a, b):
    """
    Check if the file a is newer than file b
    """
    if not (op.exists(a) and op.exists(b)):
        return False
    am = os.stat(a).st_mtime
    bm = os.stat(b).st_mtime
    return am > bm


def need_update(a, b):
    """
    Check if file a is newer than file b and decide whether or not to update
    file b. Can generalize to two lists.
    """
    if isinstance(a, basestring):
        a = [a]
    if isinstance(b, basestring):
        b = [b]

    return any((not op.exists(x)) for x in b) or \
           any(is_newer_file(x, y) for x in a for y in b)


def get_today():
    """
    Returns the date in 2010-07-14 format
    """
    from datetime import date
    return str(date.today())


def download(url, filename=None):
    from urlparse import urlsplit

    scheme, netloc, path, query, fragment = urlsplit(url)
    filename = filename or op.basename(path)
    filename = filename.strip()

    if not filename:
        filename = "index.html"

    if op.exists(filename):
        msg = "File `{0}` exists. Download skipped.".format(filename)
        logging.error(msg)
    else:
        from jcvi.utils.ez_setup import get_best_downloader

        downloader = get_best_downloader()
        downloader(url, filename)

    return filename


def getfilesize(filename, ratio=None):
    rawsize = op.getsize(filename)
    if not filename.endswith(".gz"):
        return rawsize

    import struct

    fo = open(filename, 'rb')
    fo.seek(-4, 2)
    r = fo.read()
    fo.close()
    size = struct.unpack('<I', r)[0]
    # This is only ISIZE, which is the UNCOMPRESSED modulo 2 ** 32
    if ratio is None:
        return size

    # Heuristic
    heuristicsize = rawsize / ratio
    while size < heuristicsize:
        size += 2 ** 32
    if size > 2 ** 32:
        logging.warn(\
            "Gzip file estimated uncompressed size: {0}.".format(size))

    return size


def debug():
    """
    Turn on the debugging
    """
    import logging

    from jcvi.apps.console import magenta, yellow

    format = yellow("%(asctime)s [%(module)s]")
    format += magenta(" %(message)s")
    logging.basicConfig(level=logging.DEBUG,
            format=format,
            datefmt="%H:%M:%S",
            )


def main():

    actions = (
        ('less', 'enhance the unix `less` command'),
        ('timestamp', 'record timestamps for all files in the current folder'),
        ('expand', 'move files in subfolders into the current folder'),
        ('touch', 'recover timestamps for files in the current folder'),
        ('mdownload', 'multiple download a list of files'),
        ('waitpid', 'wait for a PID to finish and then perform desired action'),
        ('notify', 'send an email/push notification'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def mdownload(args):
    """
    %prog mdownload links.txt

    Multiple download a list of files. Use formats.html.links() to extract the
    links file.
    """
    from jcvi.apps.grid import Jobs

    p = OptionParser(mdownload.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    linksfile, = args
    links = [(x.strip(),) for x in open(linksfile)]
    j = Jobs(download, links)
    j.run()


def expand(args):
    """
    %prog expand */*

    Move files in subfolders into the current folder.
    """
    debug()

    p = OptionParser(expand.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    seen = set()
    for a in args:
        oa = a.replace("/", "_")
        if oa in seen:
            logging.debug("Name collision `{0}`, ignored.".format(oa))
            continue

        cmd = "mv {0} {1}".format(a, oa)
        sh(cmd)
        seen.add(oa)


def fname():
    return sys._getframe().f_back.f_code.co_name


def get_times(filename):
    st = os.stat(filename)
    atime = st.st_atime
    mtime = st.st_mtime
    return (atime, mtime)


def timestamp(args):
    """
    %prog timestamp path > timestamp.info

    Record the timestamps for all files in the current folder.
    filename	atime	mtime

    This file can be used later to recover previous timestamps through touch().
    """
    p = OptionParser(timestamp.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    path, = args
    for root, dirs, files in os.walk(path):
        for f in files:
            filename = op.join(root, f)
            atime, mtime = get_times(filename)
            print filename, atime, mtime


def touch(args):
    """
    %prog touch timestamp.info

    Recover timestamps for files in the current folder.
    CAUTION: you must execute this in the same directory as timestamp().
    """
    from time import ctime

    p = OptionParser(touch.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    info, = args
    fp = open(info)
    for row in fp:
        path, atime, mtime = row.split()
        atime = float(atime)
        mtime = float(mtime)
        current_atime, current_mtime = get_times(path)

        # Check if the time has changed, with resolution up to 1 sec
        if int(atime) == int(current_atime) and \
           int(mtime) == int(current_mtime):
            continue

        times = [ctime(x) for x in (current_atime, current_mtime, atime, mtime)]
        msg = "{0} : ".format(path)
        msg += "({0}, {1}) => ({2}, {3})".format(*times)
        print >> sys.stderr, msg
        os.utime(path, (atime, mtime))


def snapshot(fp, p, fsize, counts=None):

    pos = int(p * fsize)
    print "==>> File `{0}`: {1} ({2}%)".format(fp.name, pos, int(p * 100))
    fp.seek(pos)
    fp.next()
    for i, row in enumerate(fp):
        if counts and i > counts:
            break
        try:
            sys.stdout.write(row)
        except IOError:
            break


def less(args):
    """
    %prog less filename position | less

    Enhance the unix `less` command by seeking to a file location first. This is
    useful to browse big files. Position is relative 0.00 - 1.00, or bytenumber.

    $ %prog less myfile 0.1      # Go to 10% of the current file and streaming
    $ %prog less myfile 0.1,0.2  # Stream at several positions
    $ %prog less myfile 100      # Go to certain byte number and streaming
    $ %prog less myfile 100,200  # Stream at several positions
    $ %prog less myfile all      # Generate a snapshot every 10% (10%, 20%, ..)
    """
    from jcvi.formats.base import must_open

    p = OptionParser(less.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    filename, pos = args
    fsize = getfilesize(filename)

    if pos == "all":
        pos = [x / 10. for x in range(0, 10)]
    else:
        pos = [float(x) for x in pos.split(",")]

    if pos[0] > 1:
        pos = [x / fsize for x in pos]

    if len(pos) > 1:
        counts = 20
    else:
        counts = None

    fp = must_open(filename)
    for p in pos:
        snapshot(fp, p, fsize, counts=counts)


# notification specific variables
valid_notif_methods = ["email"]
available_push_api = {"push" : ("pushover")}

def pushover(message, token, user, title="JCVI: Job Monitor", \
        priority=0, timestamp=None):
    """
    pushover.net python API

    <https://pushover.net/faq#library-python>
    """
    assert -1 <= priority <= 2, \
            "Priority should be and int() between -1 and 2"

    if timestamp == None:
        from time import time
        timestamp = int(time())

    retry, expire = (300, 3600) if priority == 2 \
            else (None, None)

    from httplib import HTTPSConnection
    from urllib import urlencode

    conn = HTTPSConnection("api.pushover.net:443")
    conn.request("POST", "/1/messages.json",
      urlencode({
          "token": token,
          "user": user,
          "message": message,
          "title": title,
          "priority": priority,
          "timestamp": timestamp,
          "retry": retry,
          "expire": expire,
      }), { "Content-type": "application/x-www-form-urlencoded" })
    conn.getresponse()


def pushnotify(subject, message, api="pushover", priority=0, timestamp=None):
    """
    Send push notifications using pre-existing APIs

    Requires a config `.ini` file in the user home area containing
    the necessary api tokens and user keys

    Default API: "pushover"
    """
    import types
    assert type(priority) is types.IntType and -1 <= priority <= 2, \
            "Priority should be and int() between -1 and 2"

    import ConfigParser

    cfgfile = op.join(op.expanduser("~"), "pushnotify.ini")
    Config = ConfigParser.ConfigParser()
    if op.exists(cfgfile):
        Config.read(cfgfile)
    else:
        sys.exit("Push notification config file `{0}`".format(cfgfile) + \
                 " does not exist!")

    if api == "pushover":
        cfg = ConfigSectionMap(Config, api)
        token, key = cfg["token"], cfg["user"]
        pushover(message, token, key, title=subject, \
                priority=priority, timestamp=timestamp)


def send_email(fromaddr, toaddr, subject, message):
    """
    Send an email message
    """
    from smtplib import SMTP

    SERVER = "localhost"
    message = "Subject: {0}\n{1}".format(subject, message)

    server = SMTP(SERVER)
    server.sendmail(fromaddr, toaddr, message)
    server.quit()


def notify(args):
    """
    %prog notify "Message to be sent"

    Send a message via email/push notification.
    Recipient email address is constructed by joining the login `username`
    and `dnsdomainname` of the server running the current process
    """
    valid_priorities = range(-1, 3, 1)
    valid_notif_methods.extend(available_push_api.keys())

    p = OptionParser(notify.__doc__)
    p.add_option("--method", default="email", choices=valid_notif_methods,
                 help="Specify the mode of notification [default: %default]" + \
                      " [choices: ({0})]".format(", ".join('"{0}"'.format(x) for x in valid_notif_methods)))
    p.add_option("--subject", default="JCVI: job monitor",
                 help="Specify the subject of the notification message")

    g1 = OptionGroup(p, "Optional parameters",
                "Note: Used when notification `method` is `push`")
    g1.add_option("--api", default="pushover", choices=available_push_api.values(),
                  help="Specify API used to send the push notification" + \
                  " [default: %default]")
    g1.add_option("--priority", default=0, type="int",
                  help="Message priority (-1 <= p <= 2) [default: %default]")
    g1.add_option("--timestamp", default=None, type="int", \
                  dest="timestamp", \
                  help="Message timestamp in unix format [default: %default]")
    p.add_option_group(g1)

    opts, args = p.parse_args(args)

    if len(args) == 0:
        logging.error("Please provide a brief message to be sent")
        sys.exit(not p.print_help())

    subject = opts.subject
    message = " ".join(args).strip()

    if opts.method == "email":
        user = getusername()
        domain = getdomainname()

        fromaddr = "notifier-donotreply@{0}".format(domain)
        toaddr = ["{0}@{1}".format(user, domain)]

        send_email(fromadd, toaddr, subject, message)
    else:
        pushnotify(subject, message, api=opts.api, priority=opts.priority, \
                   timestamp=opts.timestamp)


def is_running(pid):
    """Check whether pid exists in the current process table."""
    if pid < 0:
        return False
    import errno
    try:
        os.kill(pid, 0)
    except OSError, e:
        return e.errno == errno.EPERM
    else:
        return True


def waitpid(args):
    """
    %prog waitpid PID ::: ./command_to_run

    Given a PID, this script will wait for the PID to finish running and
    then perform a desired action (notify user and/or execute a new command)

    Specify `--grid` option to send the new process to the grid after waiting for PID
    """
    debug()
    from time import sleep
    from subprocess import check_output

    valid_notif_methods.extend(available_push_api.values())

    p = OptionParser(waitpid.__doc__)
    p.add_option("--notify", default=None, choices=valid_notif_methods,
                 help="Specify type of notification to be sent after waiting" + \
                      " [default: %default]" + \
                      " [choices: ({0})]".format(", ".join('"{0}"'.format(x) for x in valid_notif_methods)))
    p.add_option("--interval", default=120, type="int",
                 help="Specify interval at which PID should be monitored" + \
                      " [default: %default]")
    p.set_grid()
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    sep = ":::"
    cmd = None
    if sep in args:
        sepidx = args.index(sep)
        cmd = " ".join(args[sepidx + 1:]).strip()
        args = args[:sepidx]

    pid = int(" ".join(args).strip())

    status = is_running(pid)
    if status:
        get_origcmd = "ps -p {0} -o cmd h".format(pid)
        origcmd = check_output(get_origcmd.split()).strip()
        while is_running(pid):
            sleep(opts.interval)
    else:
        logging.debug("Process with PID {0} does not exist".format(pid))
        sys.exit()

    if opts.notify:
        hostname = check_output(["hostname"]).strip()
        method, api = opts.notify, None
        if opts.notify != "email":
            method, api = "push", opts.notify
        notifycmd = "--method={0}".format(method)
        if api:
            notifycmd += " --api={0}".format(api)
        notify(["[completed] {0}: `{1}`".format(gethostname(), origcmd), \
                notifycmd])

    if cmd is not None:
        bg = False if opts.grid else True
        sh(cmd, grid=opts.grid, background=bg)


if __name__ == '__main__':
    main()
