#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import pytest
import time


# Benchmark cblast vs blast module
@pytest.mark.benchmark(
    group="CBlastLine vs PyBlastLine", timer=time.time, disable_gc=True, warmup=False
)
def test_cblast(benchmark):
    from jcvi.formats.cblast import BlastLine

    @benchmark
    def result():
        # Test parser
        return BlastLine(
            "Os09g11510	Os08g13650	92.31	39	3	0	2273	2311	3237	3199	0.001	54.0"
        )

    assert result.hitlen == 39


@pytest.mark.benchmark(
    group="CBlastLine vs PyBlastLine", timer=time.time, disable_gc=True, warmup=False
)
def test_pyblast(benchmark):
    from jcvi.formats.pyblast import BlastLine

    @benchmark
    def result():
        # Test parser
        return BlastLine(
            "Os09g11510	Os08g13650	92.31	39	3	0	2273	2311	3237	3199	0.001	54.0"
        )

    assert result.hitlen == 39
