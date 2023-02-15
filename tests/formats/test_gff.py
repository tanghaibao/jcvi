#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_gff.py
# formats
#
# Created by Haibao Tang on 03/20/21
# Copyright © 2021 Haibao Tang. All rights reserved.
#

import pytest

from jcvi.formats.gff import DefaultOrderedDict, GffLine, make_attributes


@pytest.mark.parametrize(
    "attr,gff3,expected",
    [
        (
            "ID=cds00002;Parent=mRNA00002;",
            True,
            DefaultOrderedDict(ID=["cds00002"], Parent=["mRNA00002"]),
        ),
        (
            'Gene 22240.t000374; Note "Carbonic anhydrase"',
            False,
            DefaultOrderedDict(Gene=["22240.t000374"], Note=["Carbonic anhydrase"]),
        ),
    ],
)
def test_make_attributes(attr, gff3, expected):
    features = make_attributes(attr, gff3)
    assert features == expected


@pytest.mark.parametrize(
    "gff3_line,parent_key,expected",
    [
        (
            "mRNA00001\t.\tCDS\t1\t100\t.\t+\t.\tID=cds00001;Parent=mRNA00001;",
            "ID",
            "cds00001",
        ),
        (
            "mRNA00001\t.\tCDS\t1\t100\t.\t+\t.\tID=cds00001;Parent=mRNA00001;",
            "Parent",
            "mRNA00001",
        ),
        (
            "Chr1\tMaker\tmRNA\t851\t13757\t.\t+\t.\tID=A014204.m1;geneID=A014204;gene_name=A014204",
            "ID",
            "A014204.m1",
        ),
        (
            "Chr1\tMaker\tmRNA\t851\t13757\t.\t+\t.\tID=A014204.m1;geneID=A014204;gene_name=A014204",
            "geneID",
            "A014204",
        ),
        (
            "Chr1\tMaker\tmRNA\t851\t13757\t.\t+\t.\tID=A014204.m1;geneID=A014204;gene_name=A014204",
            "gene_name",
            "A014204",
        ),
        (
            "Chr1\tMaker\tmRNA\t851\t13757\t.\t+\t.\tID=A014204.m1;geneID=A014204;gene_name=A014204",
            "DOES_NOT_EXIST",
            None,
        ),
    ],
)
def test_parent_key(gff3_line, parent_key, expected):
    gff3_line = GffLine(gff3_line, parent_key=parent_key)
    assert gff3_line.parent == expected
