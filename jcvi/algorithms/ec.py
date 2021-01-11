#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This module contains methods to interface with DEAP evolutionary computation
framewor, including a Genetic Algorithm (GA) based method to solve scaffold
ordering and orientation problem.
"""

import sys
import array
import random
import logging
import multiprocessing

from deap import base, creator, tools
from deap.algorithms import varAnd
from jcvi.algorithms.lis import longest_monotonic_subseq_length
from jcvi.utils.console import printf


# This has to be in global space, otherwise runs into error "creator.Individual
# not found" when runnning on macOS. See also:
# https://github.com/DEAP/deap/issues/268
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", array.array, typecode="i", fitness=creator.FitnessMax)


def make_data(POINTS, SCF):
    seq = range(POINTS)
    scaffolds = []
    batch = POINTS // SCF
    for i in range(SCF):
        p = seq[i * batch : (i + 1) * batch]
        scaffolds.append(p)
    return scaffolds


def colinear_evaluate(tour, scaffolds):
    series = []
    for t in tour:
        series.extend(scaffolds[t])
    score, diff = longest_monotonic_subseq_length(series)
    return (score,)


def genome_mutation(candidate):
    """Return the mutants created by inversion mutation on the candidates.

    This function performs inversion or insertion. It randomly chooses two
    locations along the candidate and reverses the values within that
    slice. Insertion is done by popping one item and insert it back at random
    position.
    """
    size = len(candidate)
    prob = random.random()
    if prob > 0.5:  # Inversion
        p = random.randint(0, size - 1)
        q = random.randint(0, size - 1)
        if p > q:
            p, q = q, p
        q += 1
        s = candidate[p:q]
        x = candidate[:p] + s[::-1] + candidate[q:]
        return (creator.Individual(x),)
    else:  # Insertion
        p = random.randint(0, size - 1)
        q = random.randint(0, size - 1)
        cq = candidate.pop(q)
        candidate.insert(p, cq)
        return (candidate,)


def genome_mutation_orientation(candidate):
    size = len(candidate)
    prob = random.random()
    if prob > 0.5:  # Range flip
        p = random.randint(0, size - 1)
        q = random.randint(0, size - 1)
        if p > q:
            p, q = q, p
        q += 1
        for x in range(p, q):
            candidate[x] = -candidate[x]
    else:  # Single flip
        p = random.randint(0, size - 1)
        candidate[p] = -candidate[p]
    return (candidate,)


def GA_setup(guess):
    toolbox = base.Toolbox()

    toolbox.register("individual", creator.Individual, guess)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxPartialyMatched)
    toolbox.register("mutate", genome_mutation)
    toolbox.register("select", tools.selTournament, tournsize=3)
    return toolbox


def eaSimpleConverge(
    population,
    toolbox,
    cxpb,
    mutpb,
    ngen,
    stats=None,
    halloffame=None,
    callback=None,
    verbose=True,
):
    """This algorithm reproduce the simplest evolutionary algorithm as
    presented in chapter 7 of [Back2000]_.

    Modified to allow checking if there is no change for ngen, as a simple
    rule for convergence. Interface is similar to eaSimple(). However, in
    eaSimple, ngen is total number of iterations; in eaSimpleConverge, we
    terminate only when the best is NOT updated for ngen iterations.
    """
    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    if halloffame is not None:
        halloffame.update(population)

    record = stats.compile(population) if stats else {}

    # Begin the generational process
    gen = 1
    best = (0,)
    while True:
        # Select the next generation individuals
        offspring = toolbox.select(population, len(population))

        # Vary the pool of individuals
        offspring = varAnd(offspring, toolbox, cxpb, mutpb)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Update the hall of fame with the generated individuals
        if halloffame is not None:
            halloffame.update(offspring)

        if callback is not None:
            callback(halloffame[0], gen)

        # Replace the current population by the offspring
        population[:] = offspring

        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        current_best = record["max"]
        if gen % 20 == 0 and verbose:
            printf(
                "Current iteration {0}: max_score={1}".format(gen, current_best),
            )

        if current_best > best:
            best = current_best
            updated = gen

        gen += 1
        if gen - updated > ngen:
            break

    return population


def GA_run(toolbox, ngen=500, npop=100, seed=666, cpus=1, callback=None):
    logging.debug("GA setup: ngen=%d npop=%d cpus=%d seed=%d", ngen, npop, cpus, seed)
    if cpus > 1:
        pool = multiprocessing.Pool(cpus)
        toolbox.register("map", pool.map)
    random.seed(seed)
    pop = toolbox.population(n=npop)
    hof = tools.HallOfFame(1)

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("max", max)
    stats.register("min", min)

    eaSimpleConverge(
        pop, toolbox, 0.7, 0.2, ngen, stats=stats, halloffame=hof, callback=callback
    )
    tour = hof[0]
    if cpus > 1:
        pool.terminate()
    return tour, tour.fitness


if __name__ == "__main__":
    POINTS, SCF = 200, 20
    scaffolds = make_data(POINTS, SCF)

    # Demo case: scramble of the list
    guess = list(range(SCF))
    guess[5:15] = guess[5:15][::-1]
    guess[7:18] = guess[7:18][::-1]
    printf(guess)

    toolbox = GA_setup(guess)
    toolbox.register("evaluate", colinear_evaluate, scaffolds=scaffolds)
    tour, tour.fitness = GA_run(toolbox, cpus=8)
    printf(tour, tour.fitness)
