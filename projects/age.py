#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts related to age prediction model.
"""

import logging
import json
import os
import os.path as op
import sys

import pandas as pd

from jinja2 import Template
from jcvi.apps.base import OptionParser, ActionDispatcher, iglob


def main():

    actions = (
        ('compile', 'extract telomere length and ccn'),
        ('traits', 'make HTML page that reports eye and skin color'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


traits_template = '''
<html>
    <head>
        <title>ART traits</title>
        <link href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/materialize/0.100.2/css/materialize.min.css">
        <script src="https://cdnjs.cloudflare.com/ajax/libs/materialize/0.100.2/js/materialize.min.js"></script>
        <style>
        img{
            float:left;
            border-radius: 5px;
            margin-right: 10px;
            margin-left: 10px;
            width: 128;
            height: 64;
        }
        #box{
            border-radius: 50%;
            width: 48px;
            height: 48px;
        }
        </style>
    </head>
    <body class="container">
    <table class="centered bordered table-bordered">
    <thead>
    <tr>
        <th>Sample ID</th>
        <th colspan="2">Skin</th>
        <th colspan="2">Eyes</th>
    </tr>
    {% for sample in samples %}
    <tr>
        <td>{{ sample.sample_id }}</td>
        <td>
            <div id="box" style="background-color: {{ sample.skin_rgb }}"></div>
        </td>
        <td>
            <img src="{{ sample.traits['skin-color'].skin }}" />
        </td>
        <td>
            <div id="box" style="background-color: {{ sample.eye_rgb }}"></div>
        </td>
        <td>
            <img src="{{ sample.traits['eye-color'].right }}" />
            <img src="{{ sample.traits['eye-color'].left }}" />
        </td>
    </tr>
    {% endfor %}
    </thead>
    </table>
    </body>
</html>
'''

def lab2rgb(L, A, B):
    # Borrowed from:
    # <https://github.com/antimatter15/rgb-lab/blob/master/color.js>
    y = (L + 16) / 116
    x = A / 500 + y
    z = y - B / 200

    x = 0.95047 * (x * x * x if (x * x * x > 0.008856) else (x - 16 / 116) / 7.787)
    y = 1.00000 * (y * y * y if (y * y * y > 0.008856) else (y - 16 / 116) / 7.787)
    z = 1.08883 * (z * z * z if (z * z * z > 0.008856) else (z - 16 / 116) / 7.787)

    r = x * 3.2406 + y * -1.5372 + z * -0.4986
    g = x * -0.9689 + y * 1.8758 + z * 0.0415
    b = x * 0.0557 + y * -0.2040 + z * 1.0570

    r = (1.055 * r ** (1 / 2.4) - 0.055) if (r > 0.0031308) else 12.92 * r
    g = (1.055 * g ** (1 / 2.4) - 0.055) if (g > 0.0031308) else 12.92 * g
    b = (1.055 * b ** (1 / 2.4) - 0.055) if (b > 0.0031308) else 12.92 * b

    return max(0, min(1, r)) * 255, max(0, min(1, g)) * 255, max(0, min(1, b)) * 255


def make_rgb(L, A, B):
    r, g, b = lab2rgb(L, A, B)
    r = int(round(r))
    g = int(round(g))
    b = int(round(b))
    return "rgb({}, {}, {})".format(r, g, b)


def traits(args):
    """
    %prog traits directory

    Make HTML page that reports eye and skin color.
    """
    p = OptionParser(traits.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    samples = []
    for folder in args:
        targets = iglob(folder, "*-traits.json")
        if not targets:
            continue
        filename = targets[0]
        js = json.load(open(filename))
        js["skin_rgb"] = make_rgb(
            js["traits"]["skin-color"]["L"],
            js["traits"]["skin-color"]["A"],
            js["traits"]["skin-color"]["B"])
        js["eye_rgb"] = make_rgb(
            js["traits"]["eye-color"]["L"],
            js["traits"]["eye-color"]["A"],
            js["traits"]["eye-color"]["B"])
        samples.append(js)

    template = Template(traits_template)
    fw = open("report.html", "w")
    print >> fw, template.render(samples=samples)
    logging.debug("Report written to `{}`".format(fw.name))
    fw.close()


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
