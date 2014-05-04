#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This module contains a Genetic Algorithm (GA) based method to solve scaffold
ordering and orientation problem.
"""

import sys
import array
import random
import multiprocessing

from deap import base, creator, tools
from jcvi.algorithms.lis import longest_monotonous_subseq_length


def make_data(POINTS, SCF):
    seq = range(POINTS)
    scaffolds = {}
    batch = POINTS / SCF
    for i in xrange(SCF):
        s = i
        scaffolds[s] = seq[i * batch: (i + 1) * batch]
    return scaffolds


def colinear_evaluate(tour, scaffolds):
    series = []
    for t in tour:
        series.extend(scaffolds[t])
    score, diff = longest_monotonous_subseq_length(series)
    return score,


def genome_mutation(candidate):
    """Return the mutants created by inversion mutation on the candidates.

    This function performs inversion mutation. It randomly chooses two
    locations along the candidate and reverses the values within that
    slice.
    """
    size = len(candidate)
    prob = random.random()
    if prob > .7:    # Inversion
        p = random.randint(0, size-1)
        q = random.randint(0, size-1)
        if p > q:
            p, q = q, p
        q += 1
        s = candidate[p:q]
        x = candidate[:p] + s[::-1] + candidate[q:]
        return creator.Individual(x),
    else:            # Insertion
        p = random.randint(0, size-1)
        q = random.randint(0, size-1)
        cq = candidate.pop(q)
        candidate.insert(p, cq)
        return candidate,


def GA_setup(scaffolds, guess):
    ss = scaffolds.keys()
    SCF = len(scaffolds)

    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", array.array, typecode='i', fitness=creator.FitnessMax)

    toolbox = base.Toolbox()

    if guess:
        toolbox.register("individual", creator.Individual, guess)
    else:
        # Attribute generator
        toolbox.register("indices", random.sample, ss, SCF)
        # Structure initializers
        toolbox.register("individual", tools.initIterate, creator.Individual,
                                       toolbox.indices)

    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxPartialyMatched)
    toolbox.register("mutate", genome_mutation)
    toolbox.register("select", tools.selTournament, tournsize=3)
    toolbox.register("evaluate", colinear_evaluate, scaffolds=scaffolds)
    return toolbox


def eaSimpleConverge(population, toolbox, cxpb, mutpb, ngen, stats=None,
             halloffame=None, verbose=True):
    """This algorithm reproduce the simplest evolutionary algorithm as
    presented in chapter 7 of [Back2000]_.

    Modified by Haibao to allow checking if there is no change for ngen, as a
    rule for convergence.

    :param population: A list of individuals.
    :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                    operators.
    :param cxpb: The probability of mating two individuals.
    :param mutpb: The probability of mutating an individual.
    :param ngen: The number of generation.
    :param stats: A :class:`~deap.tools.Statistics` object that is updated
                  inplace, optional.
    :param halloffame: A :class:`~deap.tools.HallOfFame` object that will
                       contain the best individuals, optional.
    :param verbose: Whether or not to log the statistics.
    :returns: The final population.
    """
    from deap.algorithms import varAnd

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
    best = 0
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

        # Replace the current population by the offspring
        population[:] = offspring

        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        current_best = record['max']
        if gen % 20 == 0 and verbose:
            print >> sys.stderr, "Current iteration {0}: max_score={1}".\
                            format(gen, current_best)

        if current_best > best:
            best = current_best
            updated = gen

        gen += 1
        if gen - updated > ngen:
            break

    return population


def colinearsort(toolbox, scaffolds):
    pool = multiprocessing.Pool()
    toolbox.register("map", pool.map)
    #random.seed(666)
    pop = toolbox.population(n=100)
    hof = tools.HallOfFame(1)

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("max", max)
    stats.register("min", min)

    eaSimpleConverge(pop, toolbox, .7, .2, 1000, stats=stats,
                        halloffame=hof)
    tour = hof[0]
    return tour


if __name__ == "__main__":
    POINTS, SCF = 200, 20
    scaffolds = make_data(POINTS, SCF)

    # Demo case: scramble of the list
    guess = range(SCF)
    guess[5:15] = guess[5:15][::-1]
    guess[7:18] = guess[7:18][::-1]
    print guess

    toolbox = GA_setup(scaffolds, guess)
    tour = colinearsort(toolbox, scaffolds)
    print tour, colinear_evaluate(tour, scaffolds)
