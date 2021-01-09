#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest


@pytest.mark.parametrize(
    "input_array,expected",
    [
        (
            [0, 1, 2, 3, 4, 14, 13, 17, 16, 15, 5, 6, 7, 8, 9, 10, 11, 12, 18, 19],
            list(range(20)),
        ),
    ],
)
def test_ec(input_array, expected):
    from jcvi.algorithms.ec import (
        GA_setup,
        GA_run,
        colinear_evaluate,
        creator,
        make_data,
    )

    POINTS, SCF = 200, 20
    scaffolds = make_data(POINTS, SCF)
    toolbox = GA_setup(input_array)
    toolbox.register("evaluate", colinear_evaluate, scaffolds=scaffolds)
    tour, tour.fitness = GA_run(toolbox, cpus=8)
    print(tour, tour.fitness)

    assert list(tour) == expected
    assert tour.fitness == creator.FitnessMax((200.0,))
