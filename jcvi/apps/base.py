"""
Basic support for running library as script
"""
import errno
import os
import time
import os.path as op
import shutil
import signal
import sys
import logging
import fnmatch

from more_itertools import flatten
from http.client import HTTPSConnection
from urllib.parse import urlencode
from configparser import (
    ConfigParser,
    RawConfigParser,
    NoOptionError,
    NoSectionError,
    ParsingError,
)
from socket import gethostname
from subprocess import PIPE, call, check_call
from optparse import OptionParser as OptionP, OptionGroup, SUPPRESS_HELP

from natsort import natsorted
from rich.logging import Console, RichHandler

from jcvi import __copyright__, __version__

# http://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
nobreakbuffer = lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL)
nobreakbuffer()
os.environ["LC_ALL"] = "C"
JCVIHELP = "JCVI utility libraries {} [{}]\n".format(__version__, __copyright__)


class ActionDispatcher(object):
    """
    This class will be invoked
    a) when the base package is run via __main__, listing all MODULESs
    a) when a directory is run via __main__, listing all SCRIPTs
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
        if args[-2] == "jcvi":
            meta = "MODULE"
        elif args[-1] == "__main__":
            meta = "SCRIPT"
        else:
            meta = "ACTION"
        return meta, args

    def print_help(self):
        meta, args = self.get_meta()
        if meta == "MODULE":
            del args[0]
            args[-1] = meta
        elif meta == "SCRIPT":
            args[-1] = meta
        else:
            args[-1] += " " + meta

        help = "Usage:\n    python -m {0}\n\n\n".format(".".join(args))
        help += "Available {0}s:\n".format(meta)
        max_action_len = max(len(action) for action, ah in self.actions)
        for action, action_help in sorted(self.actions):
            action = action.rjust(max_action_len + 4)
            help += (
                " | ".join((action, action_help[0].upper() + action_help[1:])) + "\n"
            )
        help += "\n" + JCVIHELP

        sys.stderr.write(help)
        sys.exit(1)

    def dispatch(self, globals):
        from difflib import get_close_matches

        meta = "ACTION"  # function is only invoked for listing ACTIONs
        if len(sys.argv) == 1:
            self.print_help()

        action = sys.argv[1]

        if not action in self.valid_actions:
            print("[error] {0} not a valid {1}\n".format(action, meta), file=sys.stderr)
            alt = get_close_matches(action, self.valid_actions)
            print(
                "Did you mean one of these?\n\t{0}\n".format(", ".join(alt)),
                file=sys.stderr,
            )
            self.print_help()

        globals[action](sys.argv[2:])


class OptionParser(OptionP):
    def __init__(self, doc):

        OptionP.__init__(self, doc, epilog=JCVIHELP)

    def parse_args(self, args=None):
        dests = set()
        ol = []
        for g in [self] + self.option_groups:
            ol += g.option_list
        for o in ol:
            if o.dest in dests:
                continue
            self.add_help_from_choices(o)
            dests.add(o.dest)

        return OptionP.parse_args(self, args)

    def add_help_from_choices(self, o):
        if o.help == SUPPRESS_HELP:
            return

        default_tag = "%default"
        assert o.help, "Option {0} do not have help string".format(o)
        help_pf = o.help[:1].upper() + o.help[1:]
        if "[" in help_pf:
            help_pf = help_pf.rsplit("[", 1)[0]
        help_pf = help_pf.strip()

        if o.type == "choice":
            if o.default is None:
                default_tag = "guess"
            ctext = "|".join(natsorted(str(x) for x in o.choices))
            if len(ctext) > 100:
                ctext = ctext[:100] + " ... "
            choice_text = "must be one of {0}".format(ctext)
            o.help = "{0}, {1} [default: {2}]".format(help_pf, choice_text, default_tag)
        else:
            o.help = help_pf
            if o.default is None:
                default_tag = "disabled"
            if (
                o.get_opt_string() not in ("--help", "--version")
                and o.action != "store_false"
            ):
                o.help += " [default: {0}]".format(default_tag)

    def set_grid(self):
        """
        Add --grid options for command line programs
        """
        self.add_option(
            "--grid",
            dest="grid",
            default=False,
            action="store_true",
            help="Run on the grid",
        )

    def set_grid_opts(self, array=False, vcode="99999"):
        queue_choices = ("default", "fast", "medium", "himem")
        valid_pcodes = popen("qconf -sprjl", debug=False).read().strip().split("\n")
        valid_pcodes.append(vcode)

        group = OptionGroup(self, "Grid parameters")
        group.add_option(
            "-P",
            dest="pcode",
            default=vcode,
            choices=valid_pcodes,
            help="Specify accounting project code",
        )
        group.add_option(
            "-l",
            dest="queue",
            default="default",
            choices=queue_choices,
            help="Name of the queue",
        )
        group.add_option(
            "-t",
            dest="threaded",
            default=None,
            type="int",
            help="Append '-pe threaded N'",
        )
        if array:
            group.add_option(
                "-c",
                dest="concurrency",
                type="int",
                help="Append task concurrency limit '-tc N'",
            )
        group.add_option(
            "-d",
            dest="outdir",
            default=".",
            help="Specify directory to store grid output/error files",
        )
        group.add_option(
            "-N", dest="name", default=None, help="Specify descriptive name for the job"
        )
        group.add_option(
            "-H", dest="hold_jid", default=None, help="Define the job dependency list"
        )
        self.add_option_group(group)

    def set_table(self, sep=",", align=False):
        group = OptionGroup(self, "Table formatting")
        group.add_option("--sep", default=sep, help="Separator")
        if align:
            group.add_option(
                "--noalign",
                dest="align",
                default=True,
                action="store_false",
                help="Cell alignment",
            )
        else:
            group.add_option(
                "--align", default=False, action="store_true", help="Cell alignment"
            )
        self.add_option_group(group)

    def set_downloader(self, downloader=None):
        """
        Add --downloader options for given command line program.
        """
        from jcvi.utils.ez_setup import ALL_DOWNLOADERS

        downloader_choices = [x[0] for x in ALL_DOWNLOADERS]
        self.add_option(
            "--downloader",
            default=downloader,
            choices=downloader_choices,
            help="Use the specified downloader to retrieve resources",
        )

    def set_params(self, prog=None, params=""):
        """
        Add --params options for given command line programs
        """
        dest_prog = "to {0}".format(prog) if prog else ""
        self.add_option(
            "--params",
            dest="extra",
            default=params,
            help="Extra parameters to pass {0}".format(dest_prog)
            + " (these WILL NOT be validated)",
        )

    def set_outfile(self, outfile="stdout"):
        """
        Add --outfile options to print out to filename.
        """
        self.add_option("-o", "--outfile", default=outfile, help="Outfile name")

    def set_outdir(self, outdir="."):
        self.add_option("--outdir", default=outdir, help="Specify output directory")

    def set_email(self):
        """
        Add --email option to specify an email address
        """
        self.add_option(
            "--email",
            default=get_email_address(),
            help="Specify an email address",
        )

    def set_tmpdir(self, tmpdir=None):
        """
        Add --temporary_directory option to specify unix `sort` tmpdir
        """
        self.add_option(
            "-T", "--tmpdir", default=tmpdir, help="Use temp directory instead of $TMP"
        )

    def set_cpus(self, cpus=0):
        """
        Add --cpus options to specify how many threads to use.
        """
        from multiprocessing import cpu_count

        max_cpus = cpu_count()
        if not 0 < cpus < max_cpus:
            cpus = max_cpus
        self.add_option(
            "--cpus",
            default=cpus,
            type="int",
            help="Number of CPUs to use, 0=unlimited",
        )

    def set_db_opts(self, dbname="mta4", credentials=True):
        """
        Add db connection specific attributes
        """
        from jcvi.utils.db import valid_dbconn, get_profile

        self.add_option(
            "--db",
            default=dbname,
            dest="dbname",
            help="Specify name of database to query",
        )
        self.add_option(
            "--connector",
            default="Sybase",
            dest="dbconn",
            choices=valid_dbconn.keys(),
            help="Specify database connector",
        )
        hostname, username, password = get_profile()
        if credentials:
            self.add_option("--hostname", default=hostname, help="Specify hostname")
            self.add_option(
                "--username", default=username, help="Username to connect to database"
            )
            self.add_option(
                "--password", default=password, help="Password to connect to database"
            )
        self.add_option("--port", type="int", help="Specify port number")

    def set_aws_opts(self, store="hli-mv-data-science/htang"):
        from jcvi.utils.aws import s3ify

        store = s3ify(store)
        group = OptionGroup(self, "AWS and Docker options")
        self.add_option_group(group)
        # https://github.com/hlids/infrastructure/wiki/Docker-calling-convention
        group.add_option("--sample_id", help="Sample ID")
        group.add_option("--workflow_execution_id", help="Workflow execution ID")
        group.add_option("--input_bam_path", help="Input BAM location (s3 ok)")
        group.add_option("--output_path", default=store, help="Output s3 path")
        group.add_option("--workdir", default=os.getcwd(), help="Specify work dir")
        group.add_option(
            "--nocleanup",
            default=False,
            action="store_true",
            help="Don't clean up after done",
        )

    def set_stripnames(self, default=True):
        if default:
            self.add_option(
                "--no_strip_names",
                dest="strip_names",
                action="store_false",
                default=True,
                help="do not strip alternative splicing "
                "(e.g. At5g06540.1 -> At5g06540)",
            )
        else:
            self.add_option(
                "--strip_names",
                action="store_true",
                default=False,
                help="strip alternative splicing (e.g. At5g06540.1 -> At5g06540)",
            )

    def set_fixchrnames(self, orgn="medicago"):
        self.add_option(
            "--fixchrname",
            default=orgn,
            dest="fix_chr_name",
            help="Fix quirky chromosome names",
        )

    def set_SO_opts(self):
        verifySO_choices = ("verify", "resolve:prefix", "resolve:suffix")
        self.add_option(
            "--verifySO",
            choices=verifySO_choices,
            help="Verify validity of GFF3 feature type against the SO; "
            + "`resolve` will try to converge towards a valid SO "
            + "term by removing elements from the feature type "
            + "string by splitting at underscores. Example: "
            + "`mRNA_TE_gene` resolves to `mRNA` using 'resolve:prefix'",
        )

    def set_beds(self):
        self.add_option("--qbed", help="Path to qbed")
        self.add_option("--sbed", help="Path to sbed")

    def set_histogram(self, vmin=0, vmax=None, bins=20, xlabel="value", title=None):
        self.add_option(
            "--vmin", default=vmin, type="int", help="Minimum value, inclusive"
        )
        self.add_option(
            "--vmax", default=vmax, type="int", help="Maximum value, inclusive"
        )
        self.add_option(
            "--bins",
            default=bins,
            type="int",
            help="Number of bins to plot in the histogram",
        )
        self.add_option("--xlabel", default=xlabel, help="Label on the X-axis")
        self.add_option("--title", default=title, help="Title of the plot")

    def set_sam_options(self, extra=True, bowtie=False):
        self.add_option(
            "--sam",
            dest="bam",
            default=True,
            action="store_false",
            help="Write to SAM file instead of BAM",
        )
        self.add_option(
            "--uniq",
            default=False,
            action="store_true",
            help="Keep only uniquely mapped",
        )
        if bowtie:
            self.add_option(
                "--mapped", default=False, action="store_true", help="Keep mapped reads"
            )
        self.add_option(
            "--unmapped", default=False, action="store_true", help="Keep unmapped reads"
        )
        if extra:
            self.set_cpus()
            self.set_params()

    def set_mingap(self, default=100):
        self.add_option(
            "--mingap", default=default, type="int", help="Minimum size of gaps"
        )

    def set_align(
        self,
        pctid=None,
        hitlen=None,
        pctcov=None,
        evalue=None,
        compreh_pctid=None,
        compreh_pctcov=None,
        intron=None,
        bpsplice=None,
    ):
        if pctid is not None:
            self.add_option(
                "--pctid", default=pctid, type="int", help="Sequence percent identity"
            )
        if hitlen is not None:
            self.add_option(
                "--hitlen", default=hitlen, type="int", help="Minimum overlap length"
            )
        if pctcov is not None:
            self.add_option(
                "--pctcov",
                default=pctcov,
                type="int",
                help="Percentage coverage cutoff",
            )
        if evalue is not None:
            self.add_option(
                "--evalue", default=evalue, type="float", help="E-value cutoff"
            )
        if compreh_pctid is not None:
            self.add_option(
                "--compreh_pctid",
                default=compreh_pctid,
                type="int",
                help="Sequence percent identity cutoff used to "
                + "build PASA comprehensive transcriptome",
            )
        if compreh_pctcov is not None:
            self.add_option(
                "--compreh_pctcov",
                default=compreh_pctcov,
                type="int",
                help="Percent coverage cutoff used to "
                + "build PASA comprehensive transcriptome",
            )
        if intron is not None:
            self.add_option(
                "--intron",
                default=intron,
                type="int",
                help="Maximum intron length used for mapping",
            )
        if bpsplice is not None:
            self.add_option(
                "--bpsplice",
                default=bpsplice,
                type="int",
                help="Number of bp of perfect splice boundary",
            )

    def set_image_options(
        self,
        args=None,
        figsize="6x6",
        dpi=300,
        format="pdf",
        font="Helvetica",
        style="darkgrid",
        cmap="jet",
    ):
        """
        Add image format options for given command line programs.
        """
        from jcvi.graphics.base import ImageOptions, setup_theme

        allowed_format = (
            "emf",
            "eps",
            "pdf",
            "png",
            "ps",
            "raw",
            "rgba",
            "svg",
            "svgz",
        )
        allowed_fonts = ("Helvetica", "Palatino", "Schoolbook", "Arial")
        allowed_styles = ("darkgrid", "whitegrid", "dark", "white", "ticks")
        allowed_diverge = (
            "BrBG",
            "PiYG",
            "PRGn",
            "PuOr",
            "RdBu",
            "RdGy",
            "RdYlBu",
            "RdYlGn",
            "Spectral",
        )

        group = OptionGroup(self, "Image options")
        self.add_option_group(group)

        group.add_option(
            "--figsize", default=figsize, help="Figure size `width`x`height` in inches"
        )
        group.add_option(
            "--dpi",
            default=dpi,
            type="int",
            help="Physical dot density (dots per inch)",
        )
        group.add_option(
            "--format",
            default=format,
            choices=allowed_format,
            help="Generate image of format",
        )
        group.add_option(
            "--font", default=font, choices=allowed_fonts, help="Font name"
        )
        group.add_option(
            "--style", default=style, choices=allowed_styles, help="Axes background"
        )
        group.add_option(
            "--diverge",
            default="PiYG",
            choices=allowed_diverge,
            help="Contrasting color scheme",
        )
        group.add_option("--cmap", default=cmap, help="Use this color map")
        group.add_option(
            "--notex", default=False, action="store_true", help="Do not use tex"
        )

        if args is None:
            args = sys.argv[1:]

        opts, args = self.parse_args(args)

        assert opts.dpi > 0
        assert "x" in opts.figsize

        setup_theme(style=opts.style, font=opts.font, usetex=(not opts.notex))

        return opts, args, ImageOptions(opts)

    def set_dotplot_opts(self, theme=2):
        """Used in compara.catalog and graphics.dotplot"""
        from jcvi.graphics.base import set1

        group = OptionGroup(self, "Dot plot parameters")
        group.add_option(
            "--skipempty",
            default=False,
            action="store_true",
            help="Skip seqids that do not have matches",
        )
        group.add_option(
            "--nochpf",
            default=False,
            action="store_true",
            help="Do not change the contig name",
        )
        group.add_option(
            "--nostdpf",
            default=False,
            action="store_true",
            help="Do not standardize contig names",
        )
        group.add_option(
            "--genomenames",
            type="string",
            default=None,
            help="genome names for labeling axes in the form of qname_sname, "
            'eg. "*Vitis vinifera*_*Oryza sativa*"',
        )
        group.add_option(
            "--theme",
            choices=[str(x) for x in range(len(set1))],
            default=str(theme),
            help="Color index within the palette for contig grid boundaries. Palette contains: {}".format(
                "|".join(set1)
            ),
        )
        self.add_option_group(group)

    def set_depth(self, depth=50):
        self.add_option("--depth", default=depth, type="float", help="Desired depth")

    def set_rclip(self, rclip=0):
        self.add_option(
            "--rclip",
            default=rclip,
            type="int",
            help="Pair ID is derived from rstrip N chars",
        )

    def set_chr(self, chr=",".join([str(x) for x in range(1, 23)] + ["X", "Y", "MT"])):
        self.add_option("--chr", default=chr, help="Chromosomes to process")

    def set_ref(self, ref="/mnt/ref"):
        self.add_option("--ref", default=ref, help="Reference folder")

    def set_cutoff(self, cutoff=0):
        self.add_option(
            "--cutoff",
            default=cutoff,
            type="int",
            help="Distance to call valid links between mates",
        )

    def set_mateorientation(self, mateorientation=None):
        self.add_option(
            "--mateorientation",
            default=mateorientation,
            choices=("++", "--", "+-", "-+"),
            help="Use only certain mate orientations",
        )

    def set_mates(self, rclip=0, cutoff=0, mateorientation=None):
        self.set_rclip(rclip=rclip)
        self.set_cutoff(cutoff=cutoff)
        self.set_mateorientation(mateorientation=mateorientation)

    def set_bedpe(self):
        self.add_option(
            "--norc",
            dest="rc",
            default=True,
            action="store_false",
            help="Do not reverse complement, expect innie reads",
        )
        self.add_option(
            "--minlen", default=2000, type="int", help="Minimum insert size"
        )
        self.add_option(
            "--maxlen", default=8000, type="int", help="Maximum insert size"
        )
        self.add_option(
            "--dup",
            default=10,
            type="int",
            help="Filter duplicates with coordinates within this distance",
        )

    def set_fastq_names(self):
        self.add_option(
            "--names",
            default="*.fq,*.fastq,*.fq.gz,*.fastq.gz",
            help="File names to search, use comma to separate multiple",
        )

    def set_pairs(self):
        """
        %prog pairs <blastfile|samfile|bedfile>

        Report how many paired ends mapped, avg distance between paired ends, etc.
        Paired reads must have the same prefix, use --rclip to remove trailing
        part, e.g. /1, /2, or .f, .r, default behavior is to truncate until last
        char.
        """
        self.set_usage(self.set_pairs.__doc__)

        self.add_option(
            "--pairsfile", default=None, help="Write valid pairs to pairsfile"
        )
        self.add_option(
            "--nrows", default=200000, type="int", help="Only use the first n lines"
        )
        self.set_mates()
        self.add_option(
            "--pdf",
            default=False,
            action="store_true",
            help="Print PDF instead ASCII histogram",
        )
        self.add_option(
            "--bins", default=20, type="int", help="Number of bins in the histogram"
        )
        self.add_option(
            "--distmode",
            default="ss",
            choices=("ss", "ee"),
            help="Distance mode between paired reads, ss is outer distance, "
            "ee is inner distance",
        )

    def set_sep(self, sep="\t", help="Separator in the tabfile", multiple=False):
        if multiple:
            help += ", multiple values allowed"
        self.add_option("--sep", default=sep, help=help)

    def set_firstN(self, firstN=100000):
        self.add_option(
            "--firstN", default=firstN, type="int", help="Use only the first N reads"
        )

    def set_tag(self, tag=False, specify_tag=False):
        if not specify_tag:
            self.add_option(
                "--tag",
                default=tag,
                action="store_true",
                help="Add tag (/1, /2) to the read name",
            )
        else:
            tag_choices = ["/1", "/2"]
            self.add_option(
                "--tag",
                default=None,
                choices=tag_choices,
                help="Specify tag to be added to read name",
            )

    def set_phred(self, phred=None):
        phdchoices = ("33", "64")
        self.add_option(
            "--phred",
            default=phred,
            choices=phdchoices,
            help="Phred score offset {0} [default: guess]".format(phdchoices),
        )

    def set_size(self, size=0):
        self.add_option(
            "--size",
            default=size,
            type="int",
            help="Insert mean size, stdev assumed to be 20% around mean",
        )

    def set_trinity_opts(self):
        self.set_home("trinity")
        self.set_home("hpcgridrunner")
        self.set_cpus()
        self.set_params(prog="Trinity")
        topts = OptionGroup(self, "General Trinity options")
        self.add_option_group(topts)
        topts.add_option(
            "--max_memory",
            default="128G",
            type="str",
            help="Jellyfish memory allocation",
        )
        topts.add_option(
            "--min_contig_length",
            default=90,
            type="int",
            help="Minimum assembled contig length to report",
        )
        topts.add_option(
            "--bflyGCThreads",
            default=None,
            type="int",
            help="Threads for garbage collection",
        )
        topts.add_option(
            "--grid_conf_file",
            default="JCVI_SGE.0689.conf",
            type="str",
            help="HpcGridRunner config file for supported compute farms",
        )
        topts.add_option(
            "--cleanup",
            default=False,
            action="store_true",
            help="Force clean-up of unwanted files after Trinity run is complete",
        )
        ggopts = OptionGroup(self, "Genome-guided Trinity options")
        self.add_option_group(ggopts)
        ggopts.add_option(
            "--bam",
            default=None,
            type="str",
            help="provide coord-sorted bam file as starting point",
        )
        ggopts.add_option(
            "--max_intron",
            default=15000,
            type="int",
            help="maximum allowed intron length",
        )

    def set_pasa_opts(self, action="assemble"):
        self.set_home("pasa")
        if action == "assemble":
            self.set_home("tgi")
            self.add_option(
                "--clean",
                default=False,
                action="store_true",
                help="Clean transcripts using tgi seqclean",
            )
            self.set_align(pctid=95, pctcov=90, intron=15000, bpsplice=3)
            self.add_option(
                "--aligners",
                default="blat,gmap",
                help="Specify splice aligners to use for mapping",
            )
            self.add_option(
                "--fl_accs",
                default=None,
                type="str",
                help="File containing list of FL-cDNA accessions",
            )
            self.set_cpus()
            self.add_option(
                "--compreh",
                default=False,
                action="store_true",
                help="Run comprehensive transcriptome assembly",
            )
            self.set_align(compreh_pctid=95, compreh_pctcov=30)
            self.add_option(
                "--prefix",
                default="compreh_init_build",
                type="str",
                help="Prefix for compreh_trans output file names",
            )
        elif action == "compare":
            self.add_option(
                "--annots_gff3",
                default=None,
                type="str",
                help="Reference annotation to load and compare against",
            )
            genetic_code = [
                "universal",
                "Euplotes",
                "Tetrahymena",
                "Candida",
                "Acetabularia",
            ]
            self.add_option(
                "--genetic_code",
                default="universal",
                choices=genetic_code,
                help="Choose translation table",
            )
            self.add_option(
                "--pctovl",
                default=50,
                type="int",
                help="Minimum pct overlap between gene and FL assembly",
            )
            self.add_option(
                "--pct_coding",
                default=50,
                type="int",
                help="Minimum pct of cDNA sequence to be protein coding",
            )
            self.add_option(
                "--orf_size",
                default=0,
                type="int",
                help="Minimum size of ORF encoded protein",
            )
            self.add_option(
                "--utr_exons", default=2, type="int", help="Maximum number of UTR exons"
            )
            self.add_option(
                "--pctlen_FL",
                default=70,
                type="int",
                help="Minimum protein length for comparisons involving "
                + "FL assemblies",
            )
            self.add_option(
                "--pctlen_nonFL",
                default=70,
                type="int",
                help="Minimum protein length for comparisons involving "
                + "non-FL assemblies",
            )
            self.add_option(
                "--pctid_prot",
                default=70,
                type="int",
                help="Minimum pctid allowed for protein pairwise comparison",
            )
            self.add_option(
                "--pct_aln",
                default=70,
                type="int",
                help="Minimum pct of shorter protein length aligning to "
                + "update protein or isoform",
            )
            self.add_option(
                "--pctovl_gene",
                default=80,
                type="int",
                help="Minimum pct overlap among genome span of the ORF of "
                + "each overlapping gene to allow merging",
            )
            self.add_option(
                "--stompovl",
                default="",
                action="store_true",
                help="Ignore alignment results, only consider genome span of ORF",
            )
            self.add_option(
                "--trust_FL",
                default="",
                action="store_true",
                help="Trust FL-status of cDNA",
            )

    def set_annot_reformat_opts(self):
        self.add_option(
            "--pad0", default=6, type="int", help="Pad gene identifiers with 0"
        )
        self.add_option("--prefix", default="Medtr", help="Genome prefix")
        self.add_option(
            "--uc",
            default=False,
            action="store_true",
            help="Toggle gene identifier upper case",
        )

    def set_home(self, prog, default=None):
        tag = "--{0}_home".format(prog)
        default = (
            default
            or {
                "amos": "~/code/amos-code",
                "trinity": "~/export/trinityrnaseq-2.0.6",
                "hpcgridrunner": "~/export/hpcgridrunner-1.0.2",
                "cdhit": "~/export/cd-hit-v4.6.1-2012-08-27",
                "maker": "~/export/maker",
                "augustus": "~/export/maker/exe/augustus",
                "pasa": "~/export/PASApipeline-2.0.2",
                "gatk": "~/export",
                "gmes": "~/export/gmes",
                "gt": "~/export/genometools",
                "sspace": "~/export/SSPACE-STANDARD-3.0_linux-x86_64",
                "gapfiller": "~/export/GapFiller_v1-11_linux-x86_64",
                "pbjelly": "~/export/PBSuite_15.2.20",
                "picard": "~/export/picard-tools-1.138",
                "khmer": "~/export/khmer",
                "tassel": "/usr/local/projects/MTG4/packages/tassel",
                "tgi": "~/export/seqclean-x86_64",
                "eddyyeh": "/home/shared/scripts/eddyyeh",
                "fiona": "~/export/fiona-0.2.0-Linux-x86_64",
                "fermi": "~/export/fermi",
                "lobstr": "/mnt/software/lobSTR",
                "shapeit": "/mnt/software/shapeit",
                "impute": "/mnt/software/impute",
                "beagle": "java -jar /mnt/software/beagle.14Jan16.841.jar",
                "minimac": "/mnt/software/Minimac3/bin",
            }.get(prog, None)
        )
        if default is None:  # Last attempt at guessing the path
            try:
                default = op.dirname(which(prog))
            except:
                default = None
        else:
            default = op.expanduser(default)
        help = "Home directory for {0}".format(prog.upper())
        self.add_option(tag, default=default, help=help)

    def set_aligner(self, aligner="bowtie"):
        valid_aligners = ("bowtie", "bwa")
        self.add_option(
            "--aligner", default=aligner, choices=valid_aligners, help="Use aligner"
        )

    def set_verbose(self, help="Print detailed reports"):
        self.add_option("--verbose", default=False, action="store_true", help=help)


def ConfigSectionMap(Config, section):
    """
    Read a specific section from a ConfigParser() object and return
    a dict() of all key-value pairs in that section
    """
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


def get_abs_path(link_name):
    source = link_name
    if op.islink(source):
        source = os.readlink(source)
    else:
        source = op.basename(source)

    link_dir = op.dirname(link_name)
    source = op.normpath(op.join(link_dir, source))
    source = op.abspath(source)
    if source == link_name:
        return source
    else:
        return get_abs_path(source)


datadir = get_abs_path(op.join(op.dirname(__file__), "../utils/data"))
datafile = lambda x: op.join(datadir, x)


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
    """Get module-level docstring of Python module at filepath, e.g. 'path/to/file.py'."""
    co = compile(open(filepath).read(), filepath, "exec")
    if co.co_consts and isinstance(co.co_consts[0], str):
        docstring = co.co_consts[0]
    else:
        docstring = None
    return docstring


def dmain(mainfile, type="action"):
    cwd = op.dirname(mainfile)
    pyscripts = (
        [x for x in glob(op.join(cwd, "*", "__main__.py"))]
        if type == "module"
        else glob(op.join(cwd, "*.py"))
    )
    actions = []
    for ps in sorted(pyscripts):
        action = (
            op.basename(op.dirname(ps))
            if type == "module"
            else op.basename(ps).replace(".py", "")
        )
        if action[0] == "_":  # hidden namespace
            continue
        pd = get_module_docstring(ps)
        action_help = (
            [
                x.rstrip(":.,\n")
                for x in pd.splitlines(True)
                if len(x.strip()) > 10 and x[0] != "%"
            ][0]
            if pd
            else "no docstring found"
        )
        actions.append((action, action_help))

    a = ActionDispatcher(actions)
    a.print_help()


def backup(filename):
    bakname = filename + ".bak"
    if op.exists(filename):
        logging.debug("Backup `{0}` to `{1}`".format(filename, bakname))
        sh("mv {0} {1}".format(filename, bakname))
    return bakname


def getusername():
    from getpass import getuser

    return getuser()


def getdomainname():
    from socket import getfqdn

    return ".".join(str(x) for x in getfqdn().split(".")[1:])


def sh(
    cmd,
    grid=False,
    infile=None,
    outfile=None,
    errfile=None,
    append=False,
    background=False,
    threaded=None,
    log=True,
    grid_opts=None,
    silent=False,
    shell="/bin/bash",
    check=False,
):
    """
    simple wrapper for system calls
    """
    if not cmd:
        return 1
    if silent:
        outfile = errfile = "/dev/null"
    if grid:
        from jcvi.apps.grid import GridProcess

        pr = GridProcess(
            cmd,
            infile=infile,
            outfile=outfile,
            errfile=errfile,
            threaded=threaded,
            grid_opts=grid_opts,
        )
        pr.start()
        return pr.jobid
    else:
        if infile:
            cat = "cat"
            if infile.endswith(".gz"):
                cat = "zcat"
            cmd = "{0} {1} |".format(cat, infile) + cmd
        if outfile and outfile != "stdout":
            if outfile.endswith(".gz"):
                cmd += " | gzip"
            tag = ">"
            if append:
                tag = ">>"
            cmd += " {0}{1}".format(tag, outfile)
        if errfile:
            if errfile == outfile:
                errfile = "&1"
            cmd += " 2>{0}".format(errfile)
        if background:
            cmd += " &"

        if log:
            logging.debug(cmd)

        call_func = check_call if check else call
        return call_func(cmd, shell=True, executable=shell)


def Popen(cmd, stdin=None, stdout=PIPE, debug=False, shell="/bin/bash"):
    """
    Capture the cmd stdout output to a file handle.
    """
    from subprocess import Popen as P

    if debug:
        logging.debug(cmd)
    # See: <https://blog.nelhage.com/2010/02/a-very-subtle-bug/>
    proc = P(cmd, bufsize=1, stdin=stdin, stdout=stdout, shell=True, executable=shell)
    return proc


def is_macOS():
    """
    Check if current OS is macOS, this impacts mostly plotting code.
    """
    import platform

    return platform.system() == "Darwin"


def popen(cmd, debug=True, shell="/bin/bash"):
    return Popen(cmd, debug=debug, shell=shell).stdout


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


def glob(pathname, pattern=None):
    """
    Wraps around glob.glob(), but return a sorted list.
    """
    import glob as gl

    if pattern:
        pathname = op.join(pathname, pattern)
    return natsorted(gl.glob(pathname))


def iglob(pathname, patterns):
    """
    Allow multiple file formats. This is also recursive. For example:

    >>> iglob("apps", "*.py,*.pyc")
    """
    matches = []
    patterns = patterns.split(",") if "," in patterns else listify(patterns)
    for root, dirnames, filenames in os.walk(pathname):
        matching = []
        for pattern in patterns:
            matching.extend(fnmatch.filter(filenames, pattern))
        for filename in matching:
            matches.append(op.join(root, filename))
    return natsorted(matches)


def symlink(target, link_name):
    try:
        os.symlink(target, link_name)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(link_name)
            os.symlink(target, link_name)


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


def parse_multi_values(param):
    values = None
    if param:
        if op.isfile(param):
            values = list(set(x.strip() for x in open(param)))
        else:
            values = list(set(param.split(",")))
    return values


def listify(a):
    return a if (isinstance(a, list) or isinstance(a, tuple)) else [a]


def last_updated(a):
    """
    Check the time since file was last updated.
    """
    return time.time() - op.getmtime(a)


def need_update(a, b):
    """
    Check if file a is newer than file b and decide whether or not to update
    file b. Can generalize to two lists.
    """
    a = listify(a)
    b = listify(b)

    return (
        any((not op.exists(x)) for x in b)
        or all((os.stat(x).st_size == 0 for x in b))
        or any(is_newer_file(x, y) for x in a for y in b)
    )


def get_today():
    """
    Returns the date in 2010-07-14 format
    """
    from datetime import date

    return str(date.today())


def ls_ftp(dir):
    """List the contents of a remote FTP server path.

    Args:
        dir (URL): URL of a remote FTP server path.

    Returns:
        [str]: List of remote paths available, analogous to `ls`.
    """
    from urllib.parse import urlparse
    from ftpretty import ftpretty

    o = urlparse(dir)

    ftp = ftpretty(o.netloc, "anonymous", "anonymous@")
    return [op.basename(x) for x in ftp.list(o.path)]


def download(
    url, filename=None, debug=True, cookies=None, handle_gzip=False, downloader=None
):
    """Download URL to local

    Args:
        url (str): Link to the file on the internet.
        filename (str, optional): Local file name. Defaults to None.
        debug (bool, optional): Print debug messages. Defaults to True.
        cookies (str, optional): cookies file. Defaults to None.
        handle_gzip (bool, optional): Postprocess .gz files, either compress or
        uncompress. Defaults to False.
        downloader (str, optional): Use a given downloader. One of wget|curl|powershell|insecure.
        Defaults to None.

    Returns:
        str: Local file name.
    """
    from urllib.parse import urlsplit
    from subprocess import CalledProcessError
    from jcvi.formats.base import FileShredder

    scheme, netloc, path, query, fragment = urlsplit(url)
    basepath = op.basename(path)
    if basepath:
        url_gzipped = basepath.endswith(".gz")
        filename_gzipped = filename and filename.endswith(".gz")
        need_gunzip = url_gzipped and (not filename_gzipped)
        need_gzip = (not url_gzipped) and filename_gzipped
        if handle_gzip and (
            need_gunzip or need_gzip
        ):  # One more compress/decompress step after download
            target = basepath
        else:  # Just download
            target = filename or basepath
    else:
        need_gunzip, need_gzip = False, False
        target = filename or "index.html"

    success = False
    final_filename = filename or target
    if op.exists(final_filename):
        if debug:
            msg = "File `{}` exists. Download skipped.".format(final_filename)
            logging.error(msg)
        success = True
    else:
        from jcvi.utils.ez_setup import get_best_downloader

        downloader = get_best_downloader(downloader=downloader)
        if downloader:
            try:
                downloader(url, target, cookies=cookies)
                success = True
            except (CalledProcessError, KeyboardInterrupt) as e:
                print(e, file=sys.stderr)
        else:
            print("Cannot find a suitable downloader", file=sys.stderr)

        if success and handle_gzip:
            if need_gunzip:
                sh("gzip -dc {}".format(target), outfile=filename)
                FileShredder([target])
            elif need_gzip:
                sh("gzip -c {}".format(target), outfile=filename)
                FileShredder([target])

    if not success:
        FileShredder([target])

    return final_filename


def remove_if_exists(filename):
    """Check if a file exists and if so remove it

    Args:
        filename (str): Path to the local file.
    """
    if op.exists(filename):
        os.remove(filename)


def getfilesize(filename, ratio=None):
    rawsize = op.getsize(filename)
    if not filename.endswith(".gz"):
        return rawsize

    import struct

    fo = open(filename, "rb")
    fo.seek(-4, 2)
    r = fo.read()
    fo.close()
    size = struct.unpack("<I", r)[0]
    # This is only ISIZE, which is the UNCOMPRESSED modulo 2 ** 32
    if ratio is None:
        return size

    # Heuristic
    heuristicsize = rawsize / ratio
    while size < heuristicsize:
        size += 2 ** 32
    if size > 2 ** 32:
        logging.warning("Gzip file estimated uncompressed size: {0}.".format(size))

    return size


def debug(level=logging.DEBUG):
    """
    Turn on the debugging
    """
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(console=Console(stderr=True))],
    )


debug()


def main():

    actions = (
        ("less", "enhance the unix `less` command"),
        ("timestamp", "record timestamps for all files in the current folder"),
        ("expand", "move files in subfolders into the current folder"),
        ("touch", "recover timestamps for files in the current folder"),
        ("mdownload", "multiple download a list of files"),
        ("waitpid", "wait for a PID to finish and then perform desired action"),
        ("notify", "send an email/push notification"),
        ("mergecsv", "merge a set of tsv files"),
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

    (linksfile,) = args
    links = [(x.strip(),) for x in open(linksfile)]
    j = Jobs(download, links)
    j.run()


def expand(args):
    """
    %prog expand */*

    Move files in subfolders into the current folder. Use --symlink to create a
    link instead.
    """
    p = OptionParser(expand.__doc__)
    p.add_option(
        "--symlink", default=False, action="store_true", help="Create symbolic link"
    )
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    seen = set()
    for a in args:
        oa = a.replace("/", "_")
        if oa in seen:
            logging.debug("Name collision `{0}`, ignored.".format(oa))
            continue

        cmd = "cp -s" if opts.symlink else "mv"
        cmd += " {0} {1}".format(a, oa)
        sh(cmd)
        seen.add(oa)


def fname():
    return sys._getframe().f_back.f_code.co_name


def get_times(filename):
    st = os.stat(filename)
    atime = st.st_atime
    mtime = st.st_mtime
    return atime, mtime


def timestamp(args):
    """
    %prog timestamp path > timestamp.info

    Record the timestamps for all files in the current folder.
    filename atime mtime

    This file can be used later to recover previous timestamps through touch().
    """
    p = OptionParser(timestamp.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (path,) = args
    for root, dirs, files in os.walk(path):
        for f in files:
            filename = op.join(root, f)
            atime, mtime = get_times(filename)
            print(filename, atime, mtime)


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

    (info,) = args
    fp = open(info)
    for row in fp:
        path, atime, mtime = row.split()
        atime = float(atime)
        mtime = float(mtime)
        current_atime, current_mtime = get_times(path)

        # Check if the time has changed, with resolution up to 1 sec
        if int(atime) == int(current_atime) and int(mtime) == int(current_mtime):
            continue

        times = [ctime(x) for x in (current_atime, current_mtime, atime, mtime)]
        msg = "{0} : ".format(path)
        msg += "({0}, {1}) => ({2}, {3})".format(*times)
        print(msg, file=sys.stderr)
        os.utime(path, (atime, mtime))


def snapshot(fp, p, fsize, counts=None):

    pos = int(p * fsize)
    print("==>> File `{0}`: {1} ({2}%)".format(fp.name, pos, int(p * 100)))
    fp.seek(pos)
    next(fp)
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
        pos = [x / 10.0 for x in range(0, 10)]
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
available_push_api = {"push": ["pushover", "nma", "pushbullet"]}


def pushover(
    message, token, user, title="JCVI: Job Monitor", priority=0, timestamp=None
):
    """
    pushover.net python API

    <https://pushover.net/faq#library-python>
    """
    assert -1 <= priority <= 2, "Priority should be an int() between -1 and 2"

    if timestamp is None:
        from time import time

        timestamp = int(time())

    retry, expire = (300, 3600) if priority == 2 else (None, None)

    conn = HTTPSConnection("api.pushover.net:443")
    conn.request(
        "POST",
        "/1/messages.json",
        urlencode(
            {
                "token": token,
                "user": user,
                "message": message,
                "title": title,
                "priority": priority,
                "timestamp": timestamp,
                "retry": retry,
                "expire": expire,
            }
        ),
        {"Content-type": "application/x-www-form-urlencoded"},
    )
    conn.getresponse()


def nma(description, apikey, event="JCVI: Job Monitor", priority=0):
    """
    notifymyandroid.com API

    <http://www.notifymyandroid.com/api.jsp>
    """
    assert -2 <= priority <= 2, "Priority should be an int() between -2 and 2"

    conn = HTTPSConnection("www.notifymyandroid.com")
    conn.request(
        "POST",
        "/publicapi/notify",
        urlencode(
            {
                "apikey": apikey,
                "application": "python notify",
                "event": event,
                "description": description,
                "priority": priority,
            }
        ),
        {"Content-type": "application/x-www-form-urlencoded"},
    )
    conn.getresponse()


def pushbullet(body, apikey, device, title="JCVI: Job Monitor"):
    """
    pushbullet.com API

    <https://www.pushbullet.com/api>
    """
    import base64

    headers = {}
    auth = base64.encodestring("{0}:".format(apikey).encode("utf-8")).strip()
    headers["Authorization"] = "Basic {0}".format(auth)
    headers["Content-type"] = "application/x-www-form-urlencoded"

    conn = HTTPSConnection("api.pushbullet.com".format(apikey))
    conn.request(
        "POST",
        "/api/pushes",
        urlencode({"iden": device, "type": "note", "title": title, "body": body}),
        headers,
    )
    conn.getresponse()


def pushnotify(subject, message, api="pushover", priority=0, timestamp=None):
    """
    Send push notifications using pre-existing APIs

    Requires a config `pushnotify.ini` file in the user home area containing
    the necessary api tokens and user keys.

    Default API: "pushover"

    Config file format:
    -------------------
        [pushover]
        token: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        user: yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

        [nma]
        apikey: zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

        [pushbullet]
        apikey: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
        iden: dddddddddddddddddddddddddddddddddddd
    """
    assert (
        type(priority) is int and -1 <= priority <= 2
    ), "Priority should be and int() between -1 and 2"

    cfgfile = op.join(op.expanduser("~"), "pushnotify.ini")
    Config = ConfigParser()
    if op.exists(cfgfile):
        Config.read(cfgfile)
    else:
        sys.exit(
            "Push notification config file `{0}`".format(cfgfile) + " does not exist!"
        )

    if api == "pushover":
        cfg = ConfigSectionMap(Config, api)
        token, key = cfg["token"], cfg["user"]
        pushover(
            message, token, key, title=subject, priority=priority, timestamp=timestamp
        )
    elif api == "nma":
        cfg = ConfigSectionMap(Config, api)
        apikey = cfg["apikey"]
        nma(message, apikey, event=subject, priority=priority)
    elif api == "pushbullet":
        cfg = ConfigSectionMap(Config, api)
        apikey, iden = cfg["apikey"], cfg["iden"]
        pushbullet(message, apikey, iden, title=subject, type="note")


def send_email(fromaddr, toaddr, subject, message):
    """
    Send an email message
    """
    from smtplib import SMTP
    from email.mime.text import MIMEText

    SERVER = "localhost"
    _message = MIMEText(message)
    _message["Subject"] = subject
    _message["From"] = fromaddr
    _message["To"] = ", ".join(toaddr)

    server = SMTP(SERVER)
    server.sendmail(fromaddr, toaddr, _message.as_string())
    server.quit()


def get_email_address(whoami="user"):
    """Auto-generate the FROM and TO email address"""
    if whoami == "user":
        username = getusername()
        domain = getdomainname()

        myemail = "{0}@{1}".format(username, domain)
        return myemail
    else:
        fromaddr = "notifier-donotreply@{0}".format(getdomainname())
        return fromaddr


def is_valid_email(email):
    """
    RFC822 Email Address Regex
    --------------------------

    Originally written by Cal Henderson
    c.f. http://iamcal.com/publish/articles/php/parsing_email/

    Translated to Python by Tim Fletcher, with changes suggested by Dan Kubb.

    Licensed under a Creative Commons Attribution-ShareAlike 2.5 License
    http://creativecommons.org/licenses/by-sa/2.5/
    """
    import re

    qtext = "[^\\x0d\\x22\\x5c\\x80-\\xff]"
    dtext = "[^\\x0d\\x5b-\\x5d\\x80-\\xff]"
    atom = "[^\\x00-\\x20\\x22\\x28\\x29\\x2c\\x2e\\x3a-\\x3c\\x3e\\x40\\x5b-\\x5d\\x7f-\\xff]+"
    quoted_pair = "\\x5c[\\x00-\\x7f]"
    domain_literal = "\\x5b(?:%s|%s)*\\x5d" % (dtext, quoted_pair)
    quoted_string = "\\x22(?:%s|%s)*\\x22" % (qtext, quoted_pair)
    domain_ref = atom
    sub_domain = "(?:%s|%s)" % (domain_ref, domain_literal)
    word = "(?:%s|%s)" % (atom, quoted_string)
    domain = "%s(?:\\x2e%s)*" % (sub_domain, sub_domain)
    local_part = "%s(?:\\x2e%s)*" % (word, word)
    addr_spec = "%s\\x40%s" % (local_part, domain)

    email_address = re.compile(r"\A%s\Z" % addr_spec)
    if email_address.match(email):
        return True
    return False


def notify(args):
    """
    %prog notify "Message to be sent"

    Send a message via email/push notification.

    Email notify: Recipient email address is constructed by joining the login `username`
    and `dnsdomainname` of the server

    Push notify: Uses available API
    """
    valid_notif_methods.extend(available_push_api.keys())

    fromaddr = get_email_address(whoami="notifier")

    p = OptionParser(notify.__doc__)
    p.add_option(
        "--method",
        default="email",
        choices=valid_notif_methods,
        help="Specify the mode of notification",
    )
    p.add_option(
        "--subject",
        default="JCVI: job monitor",
        help="Specify the subject of the notification message",
    )
    p.set_email()

    g1 = OptionGroup(p, "Optional `push` parameters")
    g1.add_option(
        "--api",
        default="pushover",
        choices=list(flatten(available_push_api.values())),
        help="Specify API used to send the push notification",
    )
    g1.add_option(
        "--priority", default=0, type="int", help="Message priority (-1 <= p <= 2)"
    )
    g1.add_option(
        "--timestamp",
        default=None,
        type="int",
        dest="timestamp",
        help="Message timestamp in unix format",
    )
    p.add_option_group(g1)

    opts, args = p.parse_args(args)

    if len(args) == 0:
        logging.error("Please provide a brief message to be sent")
        sys.exit(not p.print_help())

    subject = opts.subject
    message = " ".join(args).strip()

    if opts.method == "email":
        toaddr = opts.email.split(",")  # TO address should be in a list
        for addr in toaddr:
            if not is_valid_email(addr):
                logging.debug("Email address `{0}` is not valid!".format(addr))
                sys.exit()
        send_email(fromaddr, toaddr, subject, message)
    else:
        pushnotify(
            subject,
            message,
            api=opts.api,
            priority=opts.priority,
            timestamp=opts.timestamp,
        )


def pid_exists(pid):
    """Check whether pid exists in the current process table."""
    if pid < 0:
        return False
    import errno

    try:
        os.kill(pid, 0)
    except OSError as e:
        return e.errno == errno.EPERM
    else:
        return True


class TimeoutExpired(Exception):
    pass


def _waitpid(pid, interval=None, timeout=None):
    """
    Wait for process with pid 'pid' to terminate and return its
    exit status code as an integer.

    If pid is not a children of os.getpid() (current process) just
    waits until the process disappears and return None.

    If pid does not exist at all return None immediately.

    Raise TimeoutExpired on timeout expired (if specified).

    Source: http://code.activestate.com/recipes/578022-wait-for-pid-and-check-for-pid-existance-posix
    """

    def check_timeout(delay):
        if timeout is not None:
            if time.time() >= stop_at:
                raise TimeoutExpired
        time.sleep(delay)
        return min(delay * 2, interval)

    if timeout is not None:
        waitcall = lambda: os.waitpid(pid, os.WNOHANG)
        stop_at = time.time() + timeout
    else:
        waitcall = lambda: os.waitpid(pid, 0)

    delay = 0.0001
    while 1:
        try:
            retpid, status = waitcall()
        except OSError as err:
            if err.errno == errno.EINTR:
                delay = check_timeout(delay)
                continue
            elif err.errno == errno.ECHILD:
                # This has two meanings:
                # - pid is not a child of os.getpid() in which case
                #   we keep polling until it's gone
                # - pid never existed in the first place
                # In both cases we'll eventually return None as we
                # can't determine its exit status code.
                while 1:
                    if pid_exists(pid):
                        delay = check_timeout(delay)
                    else:
                        return
            else:
                raise
        else:
            if retpid == 0:
                # WNOHANG was used, pid is still running
                delay = check_timeout(delay)
                continue

        # process exited due to a signal; return the integer of
        # that signal
        if os.WIFSIGNALED(status):
            return os.WTERMSIG(status)
        # process exited using exit(2) system call; return the
        # integer exit(2) system call has been called with
        elif os.WIFEXITED(status):
            return os.WEXITSTATUS(status)
        else:
            # should never happen
            raise RuntimeError("unknown process exit status")


def waitpid(args):
    """
    %prog waitpid PID ::: "./command_to_run param1 param2 ...."

    Given a PID, this script will wait for the PID to finish running and
    then perform a desired action (notify user and/or execute a new command)

    Specify "--notify=METHOD` to send the user a notification after waiting for PID
    Specify `--grid` option to send the new process to the grid after waiting for PID
    """
    import shlex

    valid_notif_methods.extend(list(flatten(available_push_api.values())))

    p = OptionParser(waitpid.__doc__)
    p.add_option(
        "--notify",
        default="email",
        choices=valid_notif_methods,
        help="Specify type of notification to be sent after waiting",
    )
    p.add_option(
        "--interval",
        default=120,
        type="int",
        help="Specify PID polling interval in seconds",
    )
    p.add_option("--message", help="Specify notification message")
    p.set_email()
    p.set_grid()
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    if not opts.message:
        """
        If notification message not specified by user, just get
        the name of the running command and use it as the message
        """
        from subprocess import check_output

    sep = ":::"
    cmd = None
    if sep in args:
        sepidx = args.index(sep)
        cmd = " ".join(args[sepidx + 1 :]).strip()
        args = args[:sepidx]

    pid = int(" ".join(args).strip())

    status = pid_exists(pid)
    if status:
        if opts.message:
            msg = opts.message
        else:
            get_origcmd = "ps -p {0} -o cmd h".format(pid)
            msg = check_output(shlex.split(get_origcmd)).strip()
        _waitpid(pid, interval=opts.interval)
    else:
        logging.debug("Process with PID {0} does not exist".format(pid))
        sys.exit()

    if opts.notify:
        notifycmd = ["[{0}] `{1}`".format(gethostname(), msg)]
        if opts.notify != "email":
            notifycmd.append("--method={0}".format("push"))
            notifycmd.append("--api={0}".format(opts.notify))
        else:
            notifycmd.append("--email={0}".format(opts.email))
        notify(notifycmd)

    if cmd is not None:
        bg = False if opts.grid else True
        sh(cmd, grid=opts.grid, background=bg)


def get_config(path):
    config = RawConfigParser()
    try:
        config.read(path)
    except ParsingError:
        e = sys.exc_info()[1]
        logging.error(
            "There was a problem reading or parsing "
            "your credentials file: %s" % (e.args[0],),
        )
    return config


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
    config = RawConfigParser()
    cfg = op.expanduser(cfg)
    changed = False
    if op.exists(cfg):
        config.read(cfg)

    assert name is not None, "Need a program name"

    try:
        fullpath = config.get(PATH, name)
    except NoSectionError:
        config.add_section(PATH)
        changed = True

    try:
        fullpath = config.get(PATH, name)
    except NoOptionError:
        msg = "=== Configure path for {0} ===\n".format(name, cfg)
        if url:
            msg += "URL: {0}\n".format(url)
        msg += "[Directory that contains `{0}`]: ".format(cmd)
        fullpath = input(msg).strip()
        config.set(PATH, name, fullpath)
        changed = True

    path = op.join(op.expanduser(fullpath), cmd)
    if warn == "exit":
        try:
            assert is_exe(path), "***ERROR: Cannot execute binary `{0}`. ".format(path)
        except AssertionError as e:
            sys.exit("{0!s}Please verify and rerun.".format(e))

    if changed:
        configfile = open(cfg, "w")
        config.write(configfile)
        logging.debug("Configuration written to `{0}`.".format(cfg))

    return path


def inspect(object):
    """A better dir() showing attributes and values"""
    for k in dir(object):
        try:
            details = getattr(object, k)
        except Exception as e:
            details = e

        try:
            details = str(details)
        except Exception as e:
            details = e

        print("{}: {}".format(k, details), file=sys.stderr)


def sample_N(a, N, seed=None):
    """When size of N is > size of a, random.sample() will emit an error:
    ValueError: sample larger than population

    This method handles such restrictions by repeatedly sampling when that
    happens. Guaranteed to cover all items if N is > size of a.

    Examples:
    >>> sample_N([1, 2, 3], 2, seed=666)
    [2, 3]
    >>> sample_N([1, 2, 3], 3, seed=666)
    [2, 3, 1]
    >>> sample_N([1, 2, 3], 4, seed=666)
    [2, 3, 1, 2]
    """
    import random

    random.seed(seed)

    ret = []
    while N > len(a):
        ret += random.sample(a, len(a))
        N -= len(a)

    return ret + random.sample(a, N)


if __name__ == "__main__":
    main()
