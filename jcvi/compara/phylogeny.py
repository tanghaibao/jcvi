#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# phylogeny.py
# compara
#
# Created by Haibao Tang on 05/21/20
# Copyright © 2020 Haibao Tang. All rights reserved.
#
import csv
import sys
import logging
import os.path as op

from jcvi.formats.fasta import Fasta, SeqIO
from jcvi.apps.base import ActionDispatcher, OptionParser, mkdir


def lcn(args):
    """
    %prog lcn Orthogroups/Orthogroups.tsv Orthogroup_Sequences/ lcn/
    """
    p = OptionParser(lcn.__doc__)
    p.add_option("--min-single-ratio", default=0.9, help="Single copy ratio must be > ")
    p.add_option("--max-zero-ratio", default=0, help="Zero copy ratio must be < ")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    (groups_tsv, sequence_dir, lcn_dir) = args
    selected = []
    # Read in the orthogroup definition and selected based on counts
    with open(groups_tsv) as fp:
        reader = csv.reader(fp, delimiter="\t")
        header = next(reader, None)
        species_names = header[1:]
        for row in reader:
            counts = [len(x.split(", ")) if x.strip() != "" else 0 for x in row[1:]]
            single_ratio = sum([x == 1 for x in counts]) / len(counts)
            zero_ratio = sum([x == 0 for x in counts]) / len(counts)
            if single_ratio < opts.min_single_ratio:
                continue
            if zero_ratio > opts.max_zero_ratio:
                continue
            print(row[0], single_ratio, zero_ratio, counts, file=sys.stderr)
            selected.append(row)

    logging.debug("A total of {} orthogroups selected".format(len(selected)))

    # Collect the FASTA sequences now
    mkdir(lcn_dir)
    for row in selected:
        orthogroup = row[0]
        orthogroup_fasta = "{}.fa".format(orthogroup)
        input_fasta = op.join(sequence_dir, orthogroup_fasta)
        fasta = Fasta(input_fasta)
        selected_seqs = []
        for gene_names, species_name in zip(row[1:], species_names):
            gene_names = gene_names.split(", ")
            if len(gene_names) == 1:
                selected, = gene_names
            else:
                max_length, selected = max((len(fasta[x]), x) for x in gene_names)
            selected_seq = fasta[selected]
            # Set gene name to species name so we can later combine them in supermatrix
            selected_seq.id = species_name
            selected_seq.name = species_name
            selected_seq.description = ""
            selected_seqs.append(selected_seq)

        output_fasta = op.join(lcn_dir, orthogroup_fasta)
        with open(output_fasta, "w") as fw:
            SeqIO.write(selected_seqs, fw, "fasta")
        print(
            "{}: {} => {} ({})".format(
                orthogroup, len(fasta), len(selected_seqs), output_fasta
            ),
            file=sys.stderr,
        )


def main():
    actions = (("lcn", "collect low copy ortholog groups from OrthoFinder results"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
