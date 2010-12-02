"""
Codes to submit multiple jobs to JCVI grid engine
"""

import sys
from subprocess import call
from optparse import OptionParser

class GridProcess (object):

    def __init__(self, cmd):
        self.cmd = cmd

    def __str__(self):
        return self.cmd

    def start(self):
        # qsub command (the project code is specific to jcvi)
        qsub = "qsub -cwd -P 04048 "

        cmd = qsub + self.cmd
        print >>sys.stderr, cmd
        call(cmd, shell=True)


class Grid (list):
    
    def __init__(self, cmds):
        for cmd in cmds:
            self.append(GridProcess(cmd))

    def run(self):
        for pi in self:
            pi.start()


def main():
    # clean the current folder, basically concatenates all prog.e* to prog.stderr
    # and prog.o* to prog.stdout
    pass


if __name__ == '__main__':
    main()
