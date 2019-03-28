#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts related to age prediction model.
"""
from __future__ import print_function

import logging
import json
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from jinja2 import Template
from jcvi.graphics.base import panel_labels, plt, savefig
import seaborn as sns
from jcvi.apps.base import OptionParser, ActionDispatcher, iglob


def main():

    actions = (
        ('compile', 'extract telomere length and ccn'),
        ('traits', 'make HTML page that reports eye and skin color'),
        # Age paper plots
        ('qc', 'plot distributions of basic statistics of a sample'),
        ('correlation', 'plot correlation of age vs. postgenomic features'),
        ('heritability', 'plot composite on heritability estimates'),
        ('regression', 'plot chronological vs. predicted age'),
        ('ccn', 'plot several ccn plots including chr1,chrX,chrY,chrM'),
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
    print(template.render(samples=samples), file=fw)
    logging.debug("Report written to `{}`".format(fw.name))
    fw.close()


def plot_fit_line(ax, x, y):
    from numpy.polynomial.polynomial import polyfit
    t = np.arange(100)
    xy = [(a, b) for (a, b) in zip(x, y) if np.isfinite(a) and np.isfinite(b)]
    x, y = zip(*xy)
    b, m = polyfit(x, y, 1)
    print("y = {} + {} * x".format(b, m))
    ax.plot(t, b + m * t, '-', lw=3, color='k')


def composite_ccn(df, size=(12, 8)):
    """ Plot composite ccn figure
    """
    fig = plt.figure(1, size)
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    ax4 = plt.subplot2grid((2, 2), (1, 1))
    chemistry = ["V1", "V2", "V2.5", float("nan")]
    colors = sns.color_palette("Set2", 8)
    color_map = dict(zip(chemistry, colors))
    mf = df[df["hli_calc_gender"] == "Male"]

    age_label = "Chronological age (yr)"
    ax1.scatter(mf["hli_calc_age_sample_taken"], mf["ccn.chrX"],
                s=10, marker='.',
                color='lightslategray')
    ax1.set_ylim(0.8, 1.1)
    plot_fit_line(ax1, mf["hli_calc_age_sample_taken"], mf["ccn.chrX"])
    ax1.set_ylabel("ChrX copy number")
    ax1.set_title("ChrX copy number in Male")

    ax2.scatter(mf["hli_calc_age_sample_taken"], mf["ccn.chrY"],
                s=10, marker='.',
                color='lightslategray')
    plot_fit_line(ax2, mf["hli_calc_age_sample_taken"], mf["ccn.chrY"])
    ax2.set_ylim(0.8, 1.1)
    ax2.set_ylabel("ChrY copy number")
    ax2.set_title("ChrY copy number in Male")

    ax3.scatter(df["hli_calc_age_sample_taken"], df["ccn.chr1"],
                s=10, marker='.',
                color='lightslategray')
    plot_fit_line(ax3, df["hli_calc_age_sample_taken"], df["ccn.chr1"])
    ax3.set_ylim(1.8, 2.1)
    ax3.set_ylabel("Chr1 copy number")
    ax3.set_title("Chr1 copy number")

    ax4.scatter(df["hli_calc_age_sample_taken"], df["ccn.chrM"],
                s=10, marker='.',
                color='lightslategray')
    plot_fit_line(ax4, df["hli_calc_age_sample_taken"], df["ccn.chrM"])
    ax4.set_ylim(0, 400)
    ax4.set_ylabel("Mitochondria copy number")
    ax4.set_title("Mitochondria copy number")

    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], marker='.', color='w', label=chem,
                          markerfacecolor=color) \
                        for (chem, color) in zip(chemistry, colors)[:3]]
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlabel(age_label)

    plt.tight_layout()
    root = fig.add_axes((0, 0, 1, 1))
    labels = ((.02, .98, "A"),
              (.52, .98, "B"),
              (.02, .5, "C"),
              (.52, .5, "D"))
    panel_labels(root, labels)
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()


def ccn(args):
    """
    %prog ccn combined.tsv

    Plot several ccn plots including chr1,chrX,chrY,chrM
    """
    p = OptionParser(ccn.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="12x8")

    if len(args) != 1:
        sys.exit(not p.print_help())

    tsvfile, = args
    df = pd.read_csv(tsvfile, sep="\t")
    composite_ccn(df, size=(iopts.w, iopts.h))
    outfile = tsvfile.rsplit(".", 1)[0] + ".ccn.pdf"
    savefig(outfile)


def regression(args):
    """
    %prog regression postgenomic-s.tsv

    Plot chronological vs. predicted age.
    """
    p = OptionParser(regression.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x8")

    if len(args) != 1:
        sys.exit(not p.print_help())

    tsvfile, = args
    df = pd.read_csv(tsvfile, sep="\t")
    chrono = "Chronological age (yr)"
    pred = "Predicted age (yr)"
    resdf = pd.DataFrame({chrono: df["hli_calc_age_sample_taken"], pred: df["Predicted Age"]})
    g = sns.jointplot(chrono, pred, resdf, joint_kws={"s": 6},
                      xlim=(0, 100), ylim=(0, 80))
    g.fig.set_figwidth(iopts.w)
    g.fig.set_figheight(iopts.h)
    outfile = tsvfile.rsplit(".", 1)[0] + ".regression.pdf"
    savefig(outfile)


def composite_correlation(df, size=(12, 8)):
    """ Plot composite correlation figure
    """
    fig = plt.figure(1, size)
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    ax4 = plt.subplot2grid((2, 2), (1, 1))
    chemistry = ["V1", "V2", "V2.5", float("nan")]
    colors = sns.color_palette("Set2", 8)
    color_map = dict(zip(chemistry, colors))

    age_label = "Chronological age (yr)"
    ax1.scatter(df["hli_calc_age_sample_taken"], df["teloLength"],
                s=10, marker='.',
                color=df["Chemistry"].map(color_map))
    ax1.set_ylim(0, 15)
    ax1.set_ylabel("Telomere length (Kb)")

    ax2.scatter(df["hli_calc_age_sample_taken"], df["ccn.chrX"],
                s=10, marker='.',
                color=df["Chemistry"].map(color_map))
    ax2.set_ylim(1.8, 2.1)
    ax2.set_ylabel("ChrX copy number")

    ax4.scatter(df["hli_calc_age_sample_taken"], df["ccn.chrY"],
                s=10, marker='.',
                color=df["Chemistry"].map(color_map))
    ax4.set_ylim(0.8, 1.1)
    ax4.set_ylabel("ChrY copy number")

    ax3.scatter(df["hli_calc_age_sample_taken"], df["TRA.PPM"],
                s=10, marker='.',
                color=df["Chemistry"].map(color_map))
    ax3.set_ylim(0, 250)
    ax3.set_ylabel("$TCR-\\alpha$ deletions (count per million reads)")

    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], marker='.', color='w', label=chem,
                          markerfacecolor=color, markersize=16) \
                        for (chem, color) in zip(chemistry, colors)[:3]]
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlabel(age_label)
        ax.legend(handles=legend_elements, loc="upper right")

    plt.tight_layout()
    root = fig.add_axes((0, 0, 1, 1))
    labels = ((.02, .98, "A"),
              (.52, .98, "B"),
              (.02, .5, "C"),
              (.52, .5, "D"))
    panel_labels(root, labels)
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()


def correlation(args):
    """
    %prog correlation postgenomic-s.tsv

    Plot correlation of age vs. postgenomic features.
    """
    p = OptionParser(correlation.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="12x8")

    if len(args) != 1:
        sys.exit(not p.print_help())

    tsvfile, = args
    df = pd.read_csv(tsvfile, sep="\t")
    composite_correlation(df, size=(iopts.w, iopts.h))
    outfile = tsvfile.rsplit(".", 1)[0] + ".correlation.pdf"
    savefig(outfile)


def composite_qc(df_orig, size=(16, 12)):
    """ Plot composite QC figures
    """
    df = df_orig.rename(columns={"hli_calc_age_sample_taken": "Age",
                       "hli_calc_gender": "Gender",
                       "eth7_max": "Ethnicity",
                       "MeanCoverage": "Mean coverage",
                       "Chemistry": "Sequencing chemistry",
                       "Release Client": "Cohort",

                      })

    fig = plt.figure(1, size)
    ax1 = plt.subplot2grid((2, 7), (0, 0), rowspan=1, colspan=2)
    ax2 = plt.subplot2grid((2, 7), (0, 2), rowspan=1, colspan=2)
    ax3 = plt.subplot2grid((2, 7), (0, 4), rowspan=1, colspan=3)
    ax4 = plt.subplot2grid((2, 7), (1, 0), rowspan=1, colspan=2)
    ax5 = plt.subplot2grid((2, 7), (1, 2), rowspan=1, colspan=2)
    ax6 = plt.subplot2grid((2, 7), (1, 4), rowspan=1, colspan=3)

    sns.distplot(df["Age"].dropna(), kde=False, ax=ax1)
    sns.countplot(x="Gender", data=df, ax=ax2)
    sns.countplot(x="Ethnicity", data=df, ax=ax3,
                    order = df['Ethnicity'].value_counts().index)
    sns.distplot(df["Mean coverage"].dropna(), kde=False, ax=ax4)
    ax4.set_xlim(0, 100)
    sns.countplot(x="Sequencing chemistry", data=df, ax=ax5)
    sns.countplot(x="Cohort", data=df, ax=ax6,
                    order = df['Cohort'].value_counts().index)
    # Anonymize the cohorts
    cohorts = ax6.get_xticklabels()
    newCohorts = []
    for i, c in enumerate(cohorts):
        if c.get_text() == "Spector":
            c = "TwinsUK"
        elif c.get_text() != "Health Nucleus":
            c = "C{}".format(i + 1)
        newCohorts.append(c)
    ax6.set_xticklabels(newCohorts)

    for ax in (ax6,):
        ax.set_xticklabels(ax.get_xticklabels(), ha="right", rotation=30)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
        ax.set_title(ax.get_xlabel())
        ax.set_xlabel("")

    plt.tight_layout()

    root = fig.add_axes((0, 0, 1, 1))
    labels = ((.02, .96, "A"),
              (.3, .96, "B"),
              (.6, .96, "C"),
              (.02, .52, "D"),
              (.3, .52, "E"),
              (.6, .52, "F"))
    panel_labels(root, labels)
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()


def qc(args):
    """
    %prog qc postgenomic-s.tsv

    Plot basic statistics of a given sample:
    Age, Gender, Ethnicity, Cohort, Chemistry
    """
    p = OptionParser(heritability.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x6")

    if len(args) != 1:
        sys.exit(not p.print_help())

    tsvfile, = args
    df = pd.read_csv(tsvfile, sep="\t")
    composite_qc(df, size=(iopts.w, iopts.h))
    outfile = tsvfile.rsplit(".", 1)[0] + ".qc.pdf"
    savefig(outfile)


def extract_trait(df, id_field, trait_field):
    traits = {}
    # Get the gender information for filtering DZ twins
    for i, row in df.iterrows():
        sample_id = str(row[id_field])
        traits[sample_id] = row[trait_field]
    return traits


def filter_same_gender(pairs, gender):
    notPresent = 0
    diffGender = 0
    # Need to screen same gender dizygotic twins
    for a, b in pairs:
        if not (a in gender and b in gender):
            notPresent += 1
            continue
        if gender[a] != gender[b]:
            diffGender += 1
            continue
        yield (a, b, gender[a])
    print(notPresent, "not found")
    print(diffGender, "different gender")


def extract_twin_values(triples, traits, gender=None):
    """Calculate the heritability of certain traits in triplets.

    Parameters
    ==========
    triples: (a, b, "Female/Male") triples. The sample IDs are then used to query
             the traits dictionary.
    traits: sample_id => value dictionary

    Returns
    =======
    tuples of size 2, that contain paired trait values of the twins
    """
    # Construct the pairs of trait values
    traitValuesAbsent = 0
    nanValues = 0
    genderSkipped = 0
    twinValues = []
    for a, b, t in triples:
        if gender is not None and t != gender:
            genderSkipped += 1
            continue
        if not (a in traits and b in traits):
            traitValuesAbsent += 1
            continue
        if np.isnan(traits[a]) or np.isnan(traits[b]):
            nanValues += 1
            continue
        twinValues.append((traits[a], traits[b]))

    print("A total of {} pairs extracted ({} absent; {} nan; {} genderSkipped)"\
        .format(len(twinValues), traitValuesAbsent, nanValues, genderSkipped))
    return twinValues


def plot_paired_values(ax, mzValues, dzValues, label=None, gender=None,
                       palette=sns.color_palette("PRGn", 10)):
    from scipy.stats import pearsonr
    mzx, mzy = zip(*mzValues)
    dzx, dzy = zip(*dzValues)
    mzline, = ax.plot(mzx, mzy, '.', color=palette[0], alpha=.75)
    dzline, = ax.plot(dzx, dzy, '.', color=palette[-1], alpha=.75)
    ax.set_xlabel(label + " in twin \#1")
    ax.set_ylabel(label + " in twin \#2")
    ax.legend((mzline, dzline), ("Monozygotic twins ($N$={}{})".format(len(mzValues), ((" " + gender) if gender else "")),
                                "Dizygotic twins ($N$={}{})".format(len(dzValues), (" " + gender) if gender else "")),
              loc="upper left")
    rho_mz, p_mz = pearsonr(mzx, mzy)
    rho_dz, p_dz = pearsonr(dzx, dzy)
    heritability = 2 * (rho_mz - rho_dz)
    ax.set_title("{} ($\\rho_{{MZ}}$={:.2f}, $\\rho_{{DZ}}$={:.2f}, $heritability$={:.2f})".\
                format(label, rho_mz, rho_dz, heritability))


def plot_abs_diff(ax, mzValues, dzValues, label=None, palette="PRGn"):
    # Let's visualize the feature differences using boxplot
    mzDelta = [abs(x - y) for (x, y) in mzValues]
    dzDelta = [abs(x - y) for (x, y) in dzValues]
    x = ["MZ twins"] * len(mzDelta) + ["DZ twins"] * len(dzDelta)
    y = mzDelta + dzDelta
    sns.boxplot(x, y, palette=palette, ax=ax)
    ax.set_ylabel("Absolute difference in {}".format(label))

def filter_low_values(data, cutoff):
    newData = [(a, b) for a, b in data if a > cutoff and b > cutoff]
    print("Removed {} outliers (<= {})"\
        .format(len(data) - len(newData), cutoff))
    return newData

def composite(df, sameGenderMZ, sameGenderDZ, size=(16, 24)):
    """Embed both absdiff figures and heritability figures.
    """
    fig = plt.figure(1, size)

    ax1a = plt.subplot2grid((6, 4), (0, 0), rowspan=2, colspan=1)
    ax2a = plt.subplot2grid((6, 4), (0, 1), rowspan=2, colspan=1)
    ax3a = plt.subplot2grid((6, 4), (0, 2), rowspan=2, colspan=1)
    ax4a = plt.subplot2grid((6, 4), (0, 3), rowspan=2, colspan=1)
    ax1b = plt.subplot2grid((6, 4), (2, 0), rowspan=2, colspan=2)
    ax2b = plt.subplot2grid((6, 4), (2, 2), rowspan=2, colspan=2)
    ax3b = plt.subplot2grid((6, 4), (4, 0), rowspan=2, colspan=2)
    ax4b = plt.subplot2grid((6, 4), (4, 2), rowspan=2, colspan=2)

    # Telomeres
    telomeres = extract_trait(df, "Sample name", "telomeres.Length")
    mzTelomeres = extract_twin_values(sameGenderMZ, telomeres)
    dzTelomeres = extract_twin_values(sameGenderDZ, telomeres)
    plot_paired_values(ax1b, mzTelomeres, dzTelomeres, label="Telomere length")
    plot_abs_diff(ax1a, mzTelomeres, dzTelomeres, label="Telomere length")

    # CCNX
    CCNX = extract_trait(df, "Sample name", "ccn.chrX")
    mzCCNX = extract_twin_values(sameGenderMZ, CCNX, gender="Female")
    dzCCNX = extract_twin_values(sameGenderDZ, CCNX, gender="Female")
    dzCCNX = filter_low_values(dzCCNX, 1.75)
    plot_paired_values(ax2b, mzCCNX, dzCCNX, gender="Female only", label="ChrX copy number")
    plot_abs_diff(ax2a, mzCCNX, dzCCNX, label="ChrX copy number")

    # CCNY
    CCNY = extract_trait(df, "Sample name", "ccn.chrY")
    mzCCNY = extract_twin_values(sameGenderMZ, CCNY, gender="Male")
    dzCCNY = extract_twin_values(sameGenderDZ, CCNY, gender="Male")
    dzCCNY = filter_low_values(dzCCNY, .75)

    plot_paired_values(ax3b, mzCCNY, dzCCNY, gender="Male only", label="ChrY copy number")
    plot_abs_diff(ax3a, mzCCNY, dzCCNY, label="ChrY copy number")

    # CCNY
    TRA = extract_trait(df, "Sample name", "TRA.PPM")
    mzTRA = extract_twin_values(sameGenderMZ, TRA)
    dzTRA = extract_twin_values(sameGenderDZ, TRA)
    plot_paired_values(ax4b, mzTRA, dzTRA, label="TCR-$\\alpha$ deletions")
    plot_abs_diff(ax4a, mzTRA, dzTRA, label="TCR-$\\alpha$ deletions")

    plt.tight_layout()

    root = fig.add_axes((0, 0, 1, 1))
    # ABCD absdiff, EFGH heritability
    labels = ((.03, .99, 'A'), (.27, .99, 'B'), (.53, .99, 'C'), (.77, .99, 'D'),
              (.03, .67, 'E'), (.53, .67, 'F'),
              (.03, .34, 'G'), (.53, .34, 'H'))
    panel_labels(root, labels)
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

def heritability(args):
    """
    %prog pg.tsv MZ-twins.csv DZ-twins.csv

    Plot composite figures ABCD on absolute difference of 4 traits,
    EFGH on heritability of 4 traits. The 4 traits are:
    telomere length, ccn.chrX, ccn.chrY, TRA.PPM
    """
    p = OptionParser(heritability.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="12x18")

    if len(args) != 3:
        sys.exit(not p.print_help())

    combined, mz, dz = args

    # Prepare twins data
    def get_pairs(filename):
        with open(filename) as fp:
            for row in fp:
                yield row.strip().split(",")

    MZ = list(get_pairs(mz))
    DZ = list(get_pairs(dz))

    print(len(MZ), "monozygotic twins")
    print(len(DZ), "dizygotic twins")

    df = pd.read_csv(combined, sep="\t", index_col=0)
    df["Sample name"] = np.array(df["Sample name"], dtype=np.str)
    gender = extract_trait(df, "Sample name", "hli_calc_gender")
    sameGenderMZ = list(filter_same_gender(MZ, gender))
    sameGenderDZ = list(filter_same_gender(DZ, gender))

    composite(df, sameGenderMZ, sameGenderDZ, size=(iopts.w, iopts.h))
    logging.getLogger().setLevel(logging.CRITICAL)
    savefig("heritability.pdf")


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
