#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Implement a few MIP solvers, based on benchmark found on <http://scip.zib.de/>
SCIP solver is ~16x faster than GLPK solver.
However, I found in rare cases it will segfault. 
Therefore the default is SCIP, the program will switch to GLPK solver for crashed cases.

The input lp_data is assumed in .lp format, see below

>>> lp_data = '''
... Maximize
...  5 x1 + 3 x2 + 2 x3
... Subject to
...  x2 + x3 <= 1
... Binary
...  x1
...  x2
...  x3
... End'''
>>> print SCIPSolver(lp_data, clean=True).results
[0, 1]
>>> print GLPKSolver(lp_data, clean=True).results
[0, 1]
"""

import os
import os.path as op
import sys
import shutil
import logging

from jcvi.apps.base import sh, debug
debug()


Work_dir = "lpsolve_work"


class AbstractMIPSolver(object):
    """
    Base class for LP solvers
    """
    def __init__(self, lp_data, work_dir=Work_dir, clean=False, verbose=False):

        self.work_dir = work_dir
        self.clean = clean
        self.verbose = verbose

        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)

        lpfile = op.join(work_dir, "data.lp") # problem instance
        logging.error("write MIP instance to `{0}`".format(lpfile))

        fw = file(lpfile, "w")
        fw.write(lp_data)
        fw.close()

        retcode, outfile = self.run(lpfile)
        if retcode < 0:
            self.results = [] 
        else:
            self.results = self.parse_output(outfile)
        
        if self.results:
            print >>sys.stderr, "optimized objective value (%d)" % self.obj_val

    def run(self, lp_data, work_dir):
        raise NotImplementedError

    def parse_output(self):
        raise NotImplementedError

    def cleanup(self):
        shutil.rmtree(self.work_dir)


class GLPKSolver(AbstractMIPSolver):
    """
    GNU Linear Programming Kit (GLPK) solver, wrapper for calling GLPSOL executable
    """
    def run(self, lpfile):

        outfile = op.join(self.work_dir, "data.lp.out") # verbose output
        listfile = op.join(self.work_dir,"data.lp.list") # simple output
        # cleanup in case something wrong happens
        for f in (outfile, listfile):
            if op.exists(f): 
                os.remove(f)

        cmd = "glpsol --cuts --fpump --lp {0} -o {1} -w {2}".format(lpfile,
                outfile, listfile)
        if not self.verbose: cmd += " >/dev/null"

        retcode = sh(cmd)

        if retcode==127:
            logging.error("You need to install program `glpsol` " + \
                          "[http://www.gnu.org/software/glpk/]")
            return -1, None

        return retcode, listfile

    def parse_output(self, listfile, clean=False):

        filtered_list = []

        fp = open(listfile)
        header = fp.readline()
        columns, rows = header.split()
        rows = int(rows)
        data = fp.readlines()
        self.obj_val = int(data[0].split()[-1])
        # the info are contained in the last several lines
        results = [int(x) for x in data[-rows:]]
        results = [i for i, x in enumerate(results) if x==1]

        fp.close()

        if self.clean: self.cleanup()

        return results


class SCIPSolver(AbstractMIPSolver):
    """
    SCIP solver, wrapper for calling SCIP executable
    """
    def run(self, lpfile):

        outfile = self.work_dir + "/data.lp.out" # verbose output
        if op.exists(outfile): 
            os.remove(outfile)

        cmd = "scip -f {0} -l {1}".format(lpfile, outfile)
        if not self.verbose:
            cmd += " >/dev/null"

        retcode = sh(cmd)

        if retcode==127:
            logging.error("You need to install program `scip` " +\
                          "[http://scip.zib.de/]")
            return -1, None

        return retcode, outfile

    def parse_output(self, outfile):

        fp = open(outfile)
        for row in fp:
            if row.startswith("objective value"): 
                obj_row = row
                break

        results = []
        for row in fp:
            #objective value:               8
            #x1                             1   (obj:5)
            #x2                             1   (obj:3)
            if row.strip()=="": break # blank line ends the section
            x = row.split()[0]
            results.append(int(x[1:])-1) # 0-based indexing

        if results:
            self.obj_val = int(obj_row.split(":")[1])

        fp.close()

        if self.clean: self.cleanup()

        return results


if __name__ == '__main__':

    import doctest
    doctest.testmod()

