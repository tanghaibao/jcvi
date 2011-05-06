#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper to call tigStore and utgcns, for debugging failed utgcns runs.

See full commands:
http://sf.net/apps/mediawiki/wgs-assembler/index.php?title=Unitig_Consensus_Failures_in_CA_6

It is expected to be executed within 5-consensus/ folder.
"""

import sys

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, sh, debug
debug()


def main():

    actions = (
        ('pull', 'pull unitig from tigStore'),
        ('test', 'test the modified unitig layout'),
        ('push', 'push the modified unitig into tigStore'),
        ('delete', 'delete specified unitig'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def pull(args):
    """
    %prog pull prefix partID unitigID

    For example,
    `%prog pull medicago 5 530` will pull the utg530 from partition 5
    The layout is written to `unitig530`
    """
    p = OptionParser(pull.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(p.print_help())

    prefix, partID, unitigID = args

    cmd = "tigStore -g ../{0}.gkpStore -t ../{0}.tigStore 1 ".format(prefix)
    cmd += "-up {0} -d layout -u {1} > unitig{0}.{1}".format(partID, unitigID)

    sh(cmd)


def test(args):
    """
    %prog test prefix partID unitigID

    For example,
    `%prog pull medicago 5 530` will test the modified `unitig530`
    """
    p = OptionParser(test.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(p.print_help())

    prefix, partID, unitigID = args

    cmd = "utgcns -g ../{0}.gkpStore -t ../{0}.tigStore 1 ".format(prefix)
    cmd += "{0} -T unitig{0}.{1} -V -V -V -v 2> unitig{0}.{1}.log".\
            format(partID, unitigID)

    sh(cmd)


def push(args):
    """
    %prog push prefix partID unitigID

    For example,
    `%prog push medicago 5 530` will push the modified `unitig530`
    and replace the one in the tigStore
    """
    p = OptionParser(push.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(p.print_help())

    prefix, partID, unitigID = args

    cmd = "tigStore -g ../{0}.gkpStore -t ../{0}.tigStore 1 ".format(prefix)
    cmd += "-up {0} -R unitig{0}.{1}".format(partID, unitigID)

    sh(cmd)


def delete(args):
    """
    %prog delete prefix partID unitigID

    For example,
    `%prog push medicago 5 530` will delete unitig 530 in partition 5
    """
    p = OptionParser(delete.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(p.print_help())

    prefix, partID, unitigID = args

    cmd = "tigStore -g ../{0}.gkpStore -t ../{0}.tigStore 1 ".format(prefix)
    cmd += "{0} -D -u {1}".format(partID, unitigID)

    sh(cmd)


if __name__ == '__main__':
    main()
