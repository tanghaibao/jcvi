"""
Pedigree file manipulation.
"""

import sys

from collections import Counter
from dataclasses import dataclass
from random import sample
from typing import Dict, Optional

import networkx as nx
import numpy as np

from ..apps.base import OptionParser, ActionDispatcher, logger, sh
from ..formats.base import BaseFile
from ..graphics.base import set3_n


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


@dataclass
class SampleInbreeding:
    """
    Store inbreeding information for a sample.
    """

    name: str
    mean_inbreeding: float
    std_inbreeding: float
    dosage: Dict[str, float]

    def __str__(self):
        return f"{self.name}\t{self.mean_inbreeding:.4f}\t{self.std_inbreeding:.4f}"


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

    def to_graph(self, inbreeding: Dict[str, SampleInbreeding]) -> nx.DiGraph:
        """
        Convert the pedigree to a graph.
        """
        G = nx.DiGraph()
        for s in self:
            dad, mom = self[s].dad, self[s].mom
            if dad:
                G.add_edge(dad, s, color="lightslategray", arrowhead="none")
            if mom:
                G.add_edge(mom, s, color="lightslategray", arrowhead="none")
        # Map colors to terminal nodes
        terminal_nodes = [s for s in self if self[s].is_terminal]
        colors = dict(zip(terminal_nodes, set3_n(len(terminal_nodes))))
        for s in self:
            inb = inbreeding[s]
            label = s
            if inb.mean_inbreeding > 0.01:
                label += f"\n(F={inb.mean_inbreeding:.2f})"
            dosage = inb.dosage
            fillcolor = [f"{colors[k]};{v:.2f}" for k, v in dosage.items()]
            fillcolor = ":".join(fillcolor)
            # Hack to make the color appear on the wedge
            if fillcolor.count(";") == 1:
                fillcolor += ":white"
            else:
                fillcolor = fillcolor.rsplit(";", 1)[0]
            G._node[s]["label"] = label
            G._node[s]["shape"] = "circle"
            G._node[s]["fixedsize"] = "true"
            G._node[s]["width"] = "0.6"
            G._node[s]["height"] = "0.6"
            G._node[s]["style"] = "wedged"
            G._node[s]["fillcolor"] = fillcolor
            G._node[s]["color"] = "none"
            G._node[s]["fontsize"] = "10"
            G._node[s]["fontname"] = "Helvetica"
        return G


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


def calculate_inbreeding(
    ped: Pedigree,
    ploidy: int,
    N: int,
) -> Dict[str, SampleInbreeding]:
    """
    Wrapper to calculate inbreeding coefficients for a sample.
    """
    logger.info("Simulating %d samples with ploidy=%d", N, ploidy)
    all_collections = []
    for _ in range(N):
        genotypes = simulate_one_iteration(ped, ploidy)
        all_collections.append(genotypes)

    results = {}
    for s in ped:
        inbreeding_coefs = [
            genotypes.inbreeding_coef(s) for genotypes in all_collections
        ]
        dosages = [genotypes.dosage(s) for genotypes in all_collections]
        dosage = sum(dosages, Counter())
        # normalize
        dosage = {k: round(v / (ploidy * N), 3) for k, v in dosage.items()}
        mean_inbreeding = float(np.mean(inbreeding_coefs))
        std_inbreeding = float(np.std(inbreeding_coefs))
        sample_inbreeding = SampleInbreeding(s, mean_inbreeding, std_inbreeding, dosage)
        results[s] = sample_inbreeding
    return results


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
    ped = Pedigree(pedfile)
    inb = calculate_inbreeding(ped, opts.ploidy, opts.N)
    print("Sample\tProportion Inbreeding\tStd dev.")
    for _, v in inb.items():
        print(v)

    G = ped.to_graph(inb)
    dotfile = f"{pedfile}.dot"
    nx.nx_agraph.write_dot(G, dotfile)
    pdf_file = dotfile + ".pdf"
    file_format = pdf_file.split(".")[-1]
    sh(f"dot -T{file_format} {dotfile} -o {pdf_file}")
    logger.info("Pedigree graph written to `%s`", pdf_file)


def main():
    actions = (("inbreeding", "calculate inbreeding coefficients"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
