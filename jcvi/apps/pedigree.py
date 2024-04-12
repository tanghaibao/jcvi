"""
Pedigree file manipulation.
"""

import sys

from dataclasses import dataclass
from typing import Optional

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

    def is_terminal(self) -> bool:
        """
        Return True if the sample is terminal.
        """
        return self.dad is None and self.mom is None


class Pedigree(BaseFile):
    """
    Read a pedigree file and store the information.
    """

    def __init__(self, pedfile: str):
        super().__init__(pedfile)
        self.samples = {}
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
                sample = Sample(name, dad, mom)
                self.samples[sample.name] = sample


def inbreeding(args):
    """
    %prog inbreeding pedfile

    Calculate inbreeding coefficients from a pedigree file.
    """
    p = OptionParser(inbreeding.__doc__)
    _, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (pedfile,) = args
    ped = Pedigree(pedfile)
    print(ped.samples)


def main():
    actions = (("inbreeding", "calculate inbreeding coefficients"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
