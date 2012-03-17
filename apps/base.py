"""
basic support for running library as script
"""

import os
import os.path as op
import sys
import shutil
import logging

from subprocess import call
from optparse import OptionParser


class ActionDispatcher (object):

    def __init__(self, actions):

        self.actions = actions
        self.valid_actions, self.action_helps = zip(*actions)

    def print_help(self):
        help = "Available actions:\n"
        for action, action_help in self.actions:
            help += "    `%s`: %s\n" % (action, action_help)

        print >>sys.stderr, help
        sys.exit(1)

    def dispatch(self, globals):
        if len(sys.argv) == 1:
            self.print_help()

        action = sys.argv[1]

        if not action in self.valid_actions:
            print >>sys.stderr, "%s not a valid action" % action
            self.print_help()

        globals[action](sys.argv[2:])


def gethostname():
    import socket
    return socket.gethostname()


def getusername():
    import getpass
    return getpass.getuser()


def set_grid(instance):
    """
    Add --grid options for command line programs
    """
    assert isinstance(instance, OptionParser)

    instance.add_option("--grid", dest="grid",
            default=False, action="store_true",
            help="run on the grid [default: %default]")


def set_params(instance):
    """
    Add --params options for given command line programs
    """
    assert isinstance(instance, OptionParser)

    instance.add_option("--params", dest="extra", default="",
            help="extra parameters to run")


def set_outfile(instance, outfile="stdout"):
    """
    Add --outfile options to print out to filename.
    """
    assert isinstance(instance, OptionParser)

    instance.add_option("-o", "--outfile", default=outfile,
            help="Outfile name [default: %default]")


def sh(cmd, grid=False, infile=None, outfile=None, errfile=None,
        background=False):
    """
    simple wrapper for system calls
    """
    if grid:
        from jcvi.apps.grid import GridProcess
        pr = GridProcess(cmd, infile=infile, outfile=outfile, errfile=errfile)
        pr.start(path=None)
        return 0  # A fake retcode
    else:
        if infile:
            cmd += " < {0} ".format(infile)
        if outfile and outfile != "stdout":
            cmd += " > {0} ".format(outfile)
        if errfile:
            cmd += " 2> {0} ".format(errfile)
        if background:
            cmd += " & "

        logging.debug(cmd)
        return call(cmd, shell=True)


def popen(cmd, debug=True):
    """
    Capture the cmd stdout output to a file handle.
    """
    from subprocess import Popen, PIPE
    if debug:
        logging.debug(cmd)
    proc = Popen(cmd, bufsize=1, stdout=PIPE, shell=True)
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
        os.mkdir(dirname)
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

    if op.exists(filename):
        msg = "File `{0}` exists. Download skipped.".format(filename)
        logging.error(msg)
    else:
        cmd = "wget {0} -O {1}".format(url, filename)
        sh(cmd)

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

    format = yellow("%(asctime)s [%(module)s::%(levelname)s]")
    format += magenta(" %(message)s")
    logging.basicConfig(level=logging.DEBUG,
            format=format,
            datefmt="%H:%M:%S",
            )


def main():

    actions = (
        ('less', 'enhance the unix `less` command'),
        ('blast', 'run blastn using query against reference'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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


def blast(args):
    """
    %prog blast ref.fasta query.fasta

    Calls blast and then filter the BLAST hits. Default is megablast.
    """
    from jcvi.apps.command import run_megablast

    task_choices = ("blastn", "blastn-short", "dc-megablast", \
                    "megablast", "vecscreen")

    p = OptionParser(blast.__doc__)
    p.add_option("--best", default=1, type="int",
            help="Only look for best N hits [default: %default]")
    p.add_option("--task", default="megablast", choices=task_choices,
            help="Task of the blastn, one of {0}".\
                 format("|".join(task_choices)) + " [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    reffasta, queryfasta = args
    q = op.basename(queryfasta).split(".")[0]
    r = op.basename(reffasta).split(".")[0]
    blastfile = "{0}.{1}.blast".format(q, r)

    run_megablast(infile=queryfasta, outfile=blastfile, db=reffasta, \
                  pctid=None, hitlen=None, best=opts.best, task=opts.task)


if __name__ == '__main__':
    main()
