#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts related to age prediction model.
"""

import json
import os
import os.path as op
import pandas as pd
import sys


from jcvi.apps.base import OptionParser, ActionDispatcher, iglob


def main():

    actions = (
        ('compile', 'extract telomere length and ccn'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def compile(args):
    """
    %prog compile directory

    Extract telomere length and ccn.
    """
    p = OptionParser(compile.__doc__)
    p.set_outfile(outfile="age.tsv")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    dfs = []
    for folder in args:
        ofolder = os.listdir(folder)

        # telomeres
        subdir = [x for x in ofolder if x.startswith("telomeres")][0]
        subdir = op.join(folder, subdir)
        filename = op.join(subdir, "tel_lengths.txt")
        df = pd.read_csv(filename, sep="\t")
        d1 = df.ix[0].to_dict()

        # ccn
        subdir = [x for x in ofolder if x.startswith("ccn")][0]
        subdir = op.join(folder, subdir)
        filename = iglob(subdir, "*.ccn.json")[0]
        js = json.load(open(filename))
        d1.update(js)
        df = pd.DataFrame(d1, index=[0])
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(opts.outfile, sep="\t", index=False)


if __name__ == '__main__':
    main()
