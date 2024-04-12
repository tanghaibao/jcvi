"""
Pedigree file manipulation.
"""

import sys

from collections import Counter
from dataclasses import dataclass
from random import sample
from typing import Optional

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

    def add(self, s: str, ploidy: int):
        """
        Add genotypes for a fixed sample (usually terminal).
        """
        self[s] = [f"{s}_{i:02d}" for i in range(ploidy)]

    def cross(self, s: str, dad: str, mom: str, ploidy: int):
        """
        Cross two samples to generate genotype for a new sample.
        """
        dad_genotype = self[dad]
        mom_genotype = self[mom]
        gamete_ploidy = ploidy // 2
        dad_gamete = sample(dad_genotype, gamete_ploidy)
        mom_gamete = sample(mom_genotype, gamete_ploidy)
        sample_genotype = sorted(dad_gamete + mom_gamete)
        self[s] = sample_genotype

    def inbreeding_coef(self, s: str) -> float:
        """
        Calculate inbreeding coefficient for a sample.

        Traditional inbreeding coefficient (F) is a measure of the probability
        that two alleles at a locus are identical by descent. This definition is
        not applicable for polyploids.

        Here we use a simpler measure of inbreeding coefficient, which is the
        proportion of alleles that are non-unique in a genotype. Or we should
        really call it "Proportion inbred".
        """
        genotype = self[s]
        ploidy = len(genotype)
        unique = len(set(genotype))
        return 1 - unique / ploidy

    def dosage(self, s: str) -> Counter:
        """
        Calculate dosage for a sample.
        """
        genotype = self[s]
        return Counter(allele.rsplit("_", 1)[0] for allele in genotype)


def simulate_one_iteration(ped: Pedigree, ploidy: int) -> GenotypeCollection:
    """
    Simulate one iteration of genotypes.
    """
    genotypes = GenotypeCollection()
    while len(genotypes) < len(ped):
        for s in ped:
            if ped[s].is_terminal:
                genotypes.add(s, ploidy=ploidy)
            else:
                dad, mom = ped[s].dad, ped[s].mom
                if dad not in genotypes or mom not in genotypes:
                    continue
                genotypes.cross(s, dad, mom, ploidy=ploidy)
    return genotypes


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
    all_collections = []
    for _ in range(N):
        genotypes = simulate_one_iteration(ped, ploidy)
        all_collections.append(genotypes)
    for s in ped:
        inbreeding_coefs = [
            genotypes.inbreeding_coef(s) for genotypes in all_collections
        ]
        dosages = [genotypes.dosage(s) for genotypes in all_collections]
        dosage = sum(dosages, Counter())
        # normalize
        dosage = {k: round(v / (ploidy * N), 3) for k, v in dosage.items()}
        mean_inbreeding = np.mean(inbreeding_coefs)
        std_inbreeding = np.std(inbreeding_coefs)
        print(f"{s}\t{mean_inbreeding:.4f}\t{std_inbreeding:.4f}\t{dosage}")


def main():
    actions = (("inbreeding", "calculate inbreeding coefficients"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
