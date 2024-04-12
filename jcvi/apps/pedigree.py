"""
Pedigree file manipulation.
"""

import sys

from collections import Counter
from dataclasses import dataclass
from random import choice, sample
from typing import Optional, Tuple

import numpy as np

from ..apps.base import OptionParser, ActionDispatcher
from ..formats.base import BaseFile


@dataclass
class Sample:
    """
    A sample in the pedigree file.
    """

    name: str
    dad: Optional[str]
    mom: Optional[str]

    @property
    def is_terminal(self) -> bool:
        """
        Return True if the sample is terminal.
        """
        return self.dad is None and self.mom is None


class Pedigree(BaseFile, dict):
    """
    Read a pedigree file and store the information.
    """

    def __init__(self, pedfile: str):
        super().__init__(pedfile)
        with open(self.filename, encoding="utf-8") as fp:
            for row in fp:
                row = row.strip()
                if row[0] == "#":  # header
                    continue
                if not row:
                    continue
                atoms = row.split()
                _, name, dad, mom = atoms[:4]
                dad = dad if dad != "0" else None
                mom = mom if mom != "0" else None
                s = Sample(name, dad, mom)
                self[s.name] = s


class GenotypeCollection(dict):
    """
    Store genotypes for each sample.
    """

    def add(self, s: str, ploidy: int, N: int):
        """
        Add genotypes for a fixed sample (usually terminal).
        """
        self[s] = [[f"{s}{i:02d}" for i in range(ploidy)] for _ in range(N)]

    def cross(self, s: str, dad: str, mom: str, ploidy: int, N: int):
        """
        Cross two samples to generate genotypes for a new sample.
        """
        dad_genotypes = self[dad]
        mom_genotypes = self[mom]
        gamete_ploidy = ploidy // 2
        sample_genotypes = []
        for _ in range(N):
            dad_genotype = choice(dad_genotypes)
            mom_genotype = choice(mom_genotypes)
            dad_gamete = sample(dad_genotype, gamete_ploidy)
            mom_gamete = sample(mom_genotype, gamete_ploidy)
            sample_genotypes.append(sorted(dad_gamete + mom_gamete))
        self[s] = sample_genotypes

    def inbreeding_coef(self, s: str) -> Tuple[float, float]:
        """
        Calculate inbreeding coefficient for a sample.
        """
        genotypes = self[s]
        results = []
        for genotype in genotypes:
            ploidy = len(genotype)
            pairs = ploidy * (ploidy - 1) // 2
            counter = Counter(genotype)
            collisions = 0
            for count in counter.values():
                collisions += count * (count - 1) // 2
            results.append(collisions / pairs)
        results = np.array(results)
        return results.mean(), results.std()


def inbreeding(args):
    """
    %prog inbreeding pedfile

    Calculate inbreeding coefficients from a pedigree file.
    """
    p = OptionParser(inbreeding.__doc__)
    p.add_option("--ploidy", default=2, type="int", help="Ploidy")
    p.add_option("--N", default=10000, type="int", help="Number of samples")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (pedfile,) = args
    ploidy = opts.ploidy
    N = opts.N
    ped = Pedigree(pedfile)
    genotypes = GenotypeCollection()
    while len(genotypes) < len(ped):
        for s in ped:
            if ped[s].is_terminal:
                genotypes.add(s, ploidy=ploidy, N=N)
            else:
                dad, mom = ped[s].dad, ped[s].mom
                if dad not in genotypes or mom not in genotypes:
                    continue
                genotypes.cross(s, dad, mom, ploidy=ploidy, N=N)
    for s in ped:
        mean, std = genotypes.inbreeding_coef(s)
        print(s, mean, std)


def main():
    actions = (("inbreeding", "calculate inbreeding coefficients"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
