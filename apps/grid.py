"""
Codes to submit multiple jobs to JCVI grid engine
"""

import sys
from subprocess import call

class Grid (object):

    def __init__(self, cmds):
        self.cmds = cmds

    def run(self):
        # qsub command (the project code is specific to jcvi)
        qsub = "qsub -cwd -P 04048 "

        for cmd in self.cmds:
            cmd = qsub + cmd
            print >>sys.stderr, cmd
            call(cmd, shell=True)
