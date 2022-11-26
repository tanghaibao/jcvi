#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
py.test script for jcvi repo

* Test cases are defined using `yml` config files (tests.yml)
* Test config files are stored in directory by the name of the script (./tests/formats/gff.py/tests.yml)
* Test directory contains valid input data files used by the test cases (./tests/formats/gff.py/inputs)
* Test directory contains reference/output files for comparison against test results (./tests/formats/gff.py/references)

This testing functionality was adapted from the implementation developed by https://github.com/CGATOxford/cgat
"""

import sys
import os.path as op
import glob
import re
import tempfile
import gzip
import time
import yaml

from importlib import import_module
from shutil import rmtree as rmdir_
from typing import Optional, Tuple

# https://stackoverflow.com/questions/16571150/how-to-capture-stdout-output-from-a-python-function-call


class Capturing(list):
    def __init__(self, stdout):
        self.fw = open(stdout, "w")

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self.fw

    def __exit__(self, *args):
        self.fw.close()
        sys.stdout = self._stdout


LOG = open("test_scripts.log", "a")


def generate_tests(metafunc, domain):
    """
    Generates a parametrized list of tests to execute
    """
    # ignore directories that don't exist as tests
    # ignore non-directories
    script_dirs = glob.glob("tests/{}/*.py".format(domain))
    script_dirs = [x for x in script_dirs if op.exists(x) and op.isdir(x)]
    script_dirs.sort()

    _test_argnames = [
        "test_name",
        "script",
        "action",
        "options",
        "arguments",
        "outputs",
        "references",
        "work_dir",
    ]
    _test_argvalues = []

    for script_dir in script_dirs:
        script_name = script_dir.replace("tests/", "")
        _script = op.abspath(script_name)
        # We need: jcvi.formats.gff
        script_name = script_name[: script_name.rfind(".py")].replace("/", ".")
        _script = "jcvi." + script_name

        yamlconf = "{0}/tests.yml".format(script_dir)
        if not op.exists(yamlconf):
            continue

        script_tests = yaml.safe_load(open(yamlconf))

        # run each test defined in SCRIPT/ACTION.py/tests.yml
        for test_name, test_params in sorted(script_tests.items()):
            _test_argvalues.append(
                (
                    test_name,
                    _script,
                    test_params["action"],
                    test_params["opts"],
                    test_params["args"],
                    test_params["outputs"],
                    test_params["references"],
                    script_dir,
                )
            )

    metafunc.parametrize(_test_argnames, _test_argvalues)


def test_script(
    test_name, script, action, options, arguments, outputs, references, work_dir
):
    """
    Runs a parametrized test with the following parameters

    :param test_name:  Unique name of the test
    :param script:     SCRIPT to invoke (e.g. `jcvi.formats.gff`)
    :param action:     SCRIPT ACTION to invoke (e.g. `format`)
    :param options:    Options to be passed (e.g. `--remove_feats=chromosome,protein`)
    :param arguments:  Arguments to be passed (e.g. `sample_input.gff`)
    :param outputs:    List of output files to collect
    :param references: List of reference files to check against output
    :param work_dir:   Directory of test data (working directory)
    :return:           Asserts test status (pass/fail)
    """
    start_time = time.time()

    tmp_dir = tempfile.mkdtemp()

    opts, args = "", ""
    if options:
        opts = _fname_resolver(options, tmp_dir=tmp_dir, work_dir=work_dir)
    if arguments:
        args = _fname_resolver(arguments, tmp_dir=tmp_dir, work_dir=work_dir)

    stdout, stderr = op.join(tmp_dir, "stdout"), op.join(tmp_dir, "stderr")

    cmd = (
        "python %(script)s %(action)s" " %(opts)s %(args)s" " > %(stdout)s"
    ) % locals()

    fail, log_msg = False, None
    module = import_module(script)
    func = getattr(module, action)
    args = " ".join((opts, args)).split()

    with Capturing(stdout) as output:
        func(args)

    if not fail:
        for output, reference in zip(outputs, references):
            output = _fname_resolver(output, tmp_dir=tmp_dir, work_dir=work_dir)

            if not op.exists(output):
                fail = True
                log_msg = "Error: output file `{0}` does not exist: {1}".format(
                    output, cmd
                )

            reference = op.join(work_dir, reference)
            if not fail and not op.exists(reference):
                fail = True
                log_msg = (
                    "Error: reference file `{0}` does not exist ({1}): {2}".format(
                        reference, tmp_dir, cmd
                    )
                )

            if not fail:
                fail, log_msg = compare_line_by_line(output, reference)

    assert not fail, log_msg

    end_time = time.time()
    LOG.write("{0}\t{1}\t{2}\n".format(script, test_name, end_time - start_time))
    LOG.flush()

    rmdir_(tmp_dir)


def compare_line_by_line(output: str, reference: str) -> Tuple[bool, Optional[str]]:
    """Compare two files line by line"""
    fail, log_msg = False, None
    for outline, refline in zip(_read(output), _read(reference)):
        if (
            outline != refline
        ):  # perform line-by-line comparison of output against reference
            fail = True
            log_msg = (
                "Error: files `{}` and `{}` are not the same\n+ {}\n- {}\n".format(
                    output, reference, outline, refline
                )
            )
            break
    return fail, log_msg


def _fname_resolver(fname, tmp_dir=None, work_dir=None):
    """
    Given a file or folder fname with placeholders (__TMP__ or __DIR__),
    replace them with values provided via named parameters

    :param tmp_dir: Temporary directory name to replace placeholder
    :param work_dir: Working directory name to replace placeholder
    :return: Returns processed file/folder name
    """
    if tmp_dir:
        fname = re.sub("__TMP__", tmp_dir, fname)
        if fname == "stdout":
            fname = op.join(tmp_dir, fname)
    if work_dir:
        fname = re.sub("__DIR__", op.abspath(work_dir), fname)

    return fname


def _read(filename):
    """
    Appropriately open file based on filetype/extension (i.e. gzip or plain-text)
    and yield lines

    :param filename: Input file (plain text or gzipped) to read
    :return: Yields lines from the input file
    """
    fp = gzip.open(filename) if filename.endswith(".gz") else open(filename)

    for line in fp:
        if not line.startswith("#"):
            yield line
