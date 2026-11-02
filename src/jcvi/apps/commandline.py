#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Command-line builders for external tools formerly wrapped by Biopython.

Recent Biopython releases no longer ship `Bio.Application` or the
`Bio.Emboss.Applications`, `Bio.Phylo.Applications` and
`Bio.Align.Applications` modules. The classes here are drop-in replacements
for the wrappers jcvi used: they accept the same keyword arguments, expose
parameters as attributes, and build the same command lines.
"""

import shlex
import subprocess
from typing import Dict, List, Optional, Tuple

from .base import logger

# (flag, keyword aliases, is_switch, equate). `equate` joins flag and value
# with "=" (EMBOSS, ClustalW style) instead of a space.
Param = Tuple[str, Tuple[str, ...], bool, bool]


def _opt(flag: str, *aliases: str, equate: bool = False) -> Param:
    return flag, aliases or (flag.lstrip("-"),), False, equate


def _switch(flag: str, *aliases: str) -> Param:
    return flag, aliases or (flag.lstrip("-"),), True, False


class Commandline:
    """
    Build and run a command line from keyword arguments.

    Subclasses set `program` (the default executable) and `params` (the
    ordered parameter table). Options with a value of None and switches that
    are falsy are omitted from the command line.
    """

    program = ""
    params: List[Param] = []

    def __init__(self, cmd: Optional[str] = None, **kwargs):
        object.__setattr__(self, "_values", {})
        object.__setattr__(self, "program", cmd or self.program)
        for name, value in kwargs.items():
            setattr(self, name, value)

    @classmethod
    def _index(cls, name: str) -> Optional[int]:
        lookup: Dict[str, int] = cls.__dict__.get("_lookup")
        if lookup is None:
            lookup = {a: i for i, p in enumerate(cls.params) for a in p[1]}
            cls._lookup = lookup
        return lookup.get(name)

    def __setattr__(self, name, value):
        idx = self._index(name)
        if idx is None:
            raise ValueError(f"{type(self).__name__} has no parameter {name!r}")
        self._values[idx] = value

    def __getattr__(self, name):
        idx = self._index(name)
        if idx is None:
            raise AttributeError(name)
        return self._values.get(idx)

    def __str__(self):
        parts = [self.program]
        for idx, (flag, _, is_switch, equate) in enumerate(self.params):
            value = self._values.get(idx)
            if is_switch:
                if value:
                    parts.append(flag)
            elif value is not None:
                value = shlex.quote(str(value))
                parts.append(f"{flag}={value}" if equate else f"{flag} {value}")
        return " ".join(parts)

    def __repr__(self):
        return f"{type(self).__name__}({str(self)!r})"

    def __call__(self) -> Tuple[str, str]:
        """
        Run the command.

        Returns:
            Tuple of (stdout, stderr). A non-zero exit status is logged as a
            warning rather than raised, so callers can check their outputs.
        """
        proc = subprocess.run(str(self), shell=True, capture_output=True, text=True)
        if proc.returncode:
            logger.warning(
                "`%s` exited with status %d: %s",
                self,
                proc.returncode,
                proc.stderr.strip()[-500:],
            )
        return proc.stdout, proc.stderr


# EMBOSS -------------------------------------------------------------------

_EMBOSS_COMMON = [
    _switch(f"-{s}")
    for s in (
        "auto",
        "stdout",
        "filter",
        "options",
        "debug",
        "verbose",
        "help",
        "warning",
        "error",
        "die",
    )
] + [_opt("-outfile", equate=True)]


def _emboss(*names: str) -> List[Param]:
    return _EMBOSS_COMMON + [_opt(f"-{n}", equate=True) for n in names]


class FSeqBootCommandline(Commandline):
    """EMBOSS PHYLIP fseqboot: bootstrap resampling of an alignment."""

    program = "fseqboot"
    params = _emboss(
        "sequence",
        "categories",
        "weights",
        "test",
        "regular",
        "fracsample",
        "rewriteformat",
        "seqtype",
        "blocksize",
        "reps",
        "justweights",
        "seed",
        "dotdiff",
    )


class FDNADistCommandline(Commandline):
    """EMBOSS PHYLIP fdnadist: nucleotide distance matrix."""

    program = "fdnadist"
    params = _emboss(
        "sequence",
        "method",
        "gamma",
        "ncategories",
        "rate",
        "categories",
        "weights",
        "gammacoefficient",
        "invarfrac",
        "ttratio",
        "freqsfrom",
        "basefreq",
        "lower",
    )


class FNeighborCommandline(Commandline):
    """EMBOSS PHYLIP fneighbor: neighbor-joining and UPGMA trees."""

    program = "fneighbor"
    params = _emboss(
        "datafile",
        "matrixtype",
        "treetype",
        "outgrno",
        "jumble",
        "seed",
        "trout",
        "outtreefile",
        "progress",
        "treeprint",
    )


class FConsenseCommandline(Commandline):
    """EMBOSS PHYLIP fconsense: majority-rule consensus tree."""

    program = "fconsense"
    params = _emboss(
        "intreefile", "method", "mlfrac", "root", "outgrno", "trout", "outtreefile"
    )


class NeedleCommandline(Commandline):
    """EMBOSS needle: Needleman-Wunsch global alignment."""

    program = "needle"
    params = (
        _emboss(
            "asequence",
            "bsequence",
            "gapopen",
            "gapextend",
            "datafile",
            "endweight",
            "endopen",
            "endextend",
        )
        + [_switch("-nobrief"), _switch("-brief")]
        + [
            _opt(f"-{n}", equate=True)
            for n in ("similarity", "snucleotide", "sprotein", "aformat")
        ]
    )


# Phylogenetics ------------------------------------------------------------


class PhymlCommandline(Commandline):
    """PhyML: maximum likelihood phylogeny."""

    program = "phyml"
    params = [
        _opt("-i", "input"),
        _opt("-d", "datatype"),
        _switch("-q", "sequential"),
        _opt("-n", "multiple"),
        _switch("-p", "pars"),
        _opt("-b", "bootstrap"),
        _opt("-m", "model"),
        _opt("-f", "frequencies"),
        _opt("-t", "ts_tv_ratio"),
        _opt("-v", "prop_invar"),
        _opt("-c", "nclasses"),
        _opt("-a", "alpha"),
        _opt("-s", "search"),
        _opt("-u", "input_tree"),
        _opt("-o", "optimize"),
        _switch("--rand_start", "rand_start"),
        _opt("--n_rand_starts", "n_rand_starts"),
        _opt("--r_seed", "r_seed"),
        _switch("--print_site_lnl", "print_site_lnl"),
        _switch("--print_trace", "print_trace"),
        _opt("--run_id", "run_id"),
        _switch("--quiet", "quiet"),
    ]


class RaxmlCommandline(Commandline):
    """RAxML: maximum likelihood phylogeny. `parsimony_seed` defaults to 10000."""

    program = "raxmlHPC"
    params = [
        _opt("-a", "weight_filename"),
        _opt("-b", "bootstrap_seed"),
        _opt("-c", "num_categories"),
        _switch("-d", "random_starting_tree"),
        _opt("-e", "epsilon"),
        _opt("-E", "exclude_filename"),
        _opt("-f", "algorithm"),
        _opt("-g", "grouping_constraint"),
        _opt("-i", "rearrangements"),
        _switch("-j", "checkpoints"),
        _switch("-k", "bootstrap_branch_lengths"),
        _opt("-l", "cluster_threshold"),
        _opt("-L", "cluster_threshold_fast"),
        _opt("-m", "model"),
        _switch("-M", "partition_branch_lengths"),
        _opt("-n", "name"),
        _opt("-o", "outgroup"),
        _opt("-q", "partition_filename"),
        _opt("-p", "parsimony_seed"),
        _opt("-P", "protein_model"),
        _opt("-r", "binary_constraint"),
        _opt("-s", "sequences"),
        _opt("-t", "starting_tree"),
        _opt("-T", "threads"),
        _opt("-u", "num_bootstrap_searches"),
        _switch("-v", "version"),
        _opt("-w", "working_dir"),
        _opt("-x", "rapid_bootstrap_seed"),
        _switch("-y", "parsimony"),
        _opt("-z", "bipartition_filename"),
        _opt("-N", "num_replicates"),
    ]

    def __init__(self, cmd: Optional[str] = None, **kwargs):
        super().__init__(cmd, **kwargs)
        if not self.parsimony_seed:
            self.parsimony_seed = 10000


# Multiple sequence alignment ----------------------------------------------

_CLUSTALW_PARAMS = """
infile profile1 profile2 options help check fullhelp align tree pim
bootstrap convert quicktree type negative outfile output outorder case
seqnos seqno_range range maxseqlen quiet stats ktuple topdiags window
pairgap score pwmatrix pwdnamatrix pwgapopen pwgapext newtree usetree matrix
dnamatrix gapopen gapext endgaps gapdist nopgap nohgap hgapresidues maxdiv
transweight iteration numiter noweights profile newtree1 newtree2 usetree1
usetree2 sequences nosecstr1 nosecstr2 secstrout helixgap strandgap loopgap
terminalgap helixendin helixendout strandendin strandendout outputtree seed
kimura tossgaps bootlabels clustering
""".split()
_CLUSTALW_SWITCHES = set("""
options help check fullhelp align tree pim convert quicktree negative quiet
endgaps nopgap nohgap hgapresidues noweights profile sequences nosecstr1
nosecstr2 kimura tossgaps
""".split())

_MUSCLE_PARAMS = """
in out diags profile in1 in2 anchorspacing center cluster1 cluster2
diaglength diagmargin distance1 distance2 gapextend gapopen hydro
hydrofactor log loga matrix diagbreak maxdiagbreak maxhours maxiters
maxtrees minbestcolscore minsmoothscore objscore refinewindow root1 root2
scorefile seqtype smoothscoreceil smoothwindow spscore sueff tree1 tree2
usetree weight1 weight2 clw clwstrict fasta html msf phyi phys phyiout
physout htmlout clwout clwstrictout msfout fastaout anchors noanchors
brenner cluster dimer group le sv sp spn quiet refine refinew core nocore
stable verbose version
""".split()
_MUSCLE_SWITCHES = set("""
diags profile clw clwstrict fasta html msf phyi phys anchors noanchors
brenner cluster dimer group le sv sp spn quiet refine refinew core nocore
stable verbose version
""".split())


class ClustalwCommandline(Commandline):
    """ClustalW multiple sequence alignment. Keywords are case-insensitive."""

    program = "clustalw"
    params = [
        (
            _switch(f"-{n}", n, n.upper())
            if n in _CLUSTALW_SWITCHES
            else _opt(f"-{n}", n, n.upper(), equate=True)
        )
        for n in _CLUSTALW_PARAMS
    ]


class MuscleCommandline(Commandline):
    """MUSCLE multiple sequence alignment. Pass `input` for the `-in` flag."""

    program = "muscle"
    params = [
        (
            _switch(f"-{n}")
            if n in _MUSCLE_SWITCHES
            else _opt(f"-{n}", *((n, "input") if n == "in" else (n,)))
        )
        for n in _MUSCLE_PARAMS
    ]
