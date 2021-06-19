#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Based on AMOS `amos-3.0.0/src/PythonModules`

AMOS message reader/writer
Contributed by Paul Harrison
"""
from __future__ import print_function

import re
import sys
import logging

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcvi.utils.cbook import percentage
from jcvi.apps.base import OptionParser, ActionDispatcher


class Message(list):
    """
    AMOS Message object

    Fields:
        type         - message type code
        fields       - dictionary of fields
        sub_messages - list of sub-messages

    str(message) converts the message back to AMOS format
    """

    def __init__(self, type):
        self.type = type
        self.contents = []

    def __str__(self):
        result = "{" + self.type + "\n"
        for key, value, multiline in self.contents:
            result += key + ":"
            if multiline:
                result += "\n{0}.\n".format(value)
            else:
                result += value + "\n"

        result += "\n".join(str(sub_message) for sub_message in self) + "}"

        return result

    def get_field(self, field):
        for key, value, multiline in self.contents:
            if key == field:
                return value
        assert ValueError("Field `{0}` cannot be found.".format(field))


_START = re.compile(r"^{([A-Z][A-Z][A-Z])\n$")
_MULTILINE_FIELD = re.compile(r"^([a-z][a-z][a-z]):\n$")
_FIELD = re.compile(r"^([a-z][a-z][a-z]):(.*)\n$")


def read_record(fp, first_line=None):
    """
    Read a record from a file of AMOS messages

    On success returns a Message object
    On end of file raises EOFError
    """

    if first_line is None:
        first_line = fp.readline()

    if not first_line:
        raise EOFError()

    match = _START.match(first_line)
    if not match:
        raise Exception("Bad start of message", first_line)

    type = match.group(1)
    message = Message(type)

    while True:
        row = fp.readline()

        match = _MULTILINE_FIELD.match(row)
        if match:
            key = match.group(1)
            val = ""
            while row:
                pos = fp.tell()
                row = fp.readline()
                if row[0] in ".":
                    break
                elif row[0] in "{}":
                    fp.seek(pos)  # put the line back
                    break
                val += row
            message.contents.append((key, val, True))
            continue

        match = _FIELD.match(row)
        if match:
            key, val = match.group(1), match.group(2)
            message.contents.append((key, val, False))
            continue

        match = _START.match(row)
        if match:
            message.append(read_record(fp, row))
            continue

        if row[0] == "}":
            break

        raise Exception("Bad line", row)

    return message


def iter_records(file):
    """Iterate over all the records in a file"""

    while True:
        try:
            yield read_record(file)
        except EOFError:
            break


def main():

    actions = (
        ("frg", "extract fasta sequences from frg"),
        ("asm", "extract fasta sequences from asm"),
        ("filter", "remove duplicate reads from frg"),
        ("count", "count each type of messages"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def filter(args):
    """
    %prog filter frgfile idsfile

    Removes the reads from frgfile that are indicated as duplicates in the
    clstrfile (generated by CD-HIT-454). `idsfile` includes a set of names to
    include in the filtered frgfile. See apps.cdhit.ids().
    """
    p = OptionParser(filter.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    frgfile, idsfile = args
    assert frgfile.endswith(".frg")

    fp = open(idsfile)
    allowed = set(x.strip() for x in fp)
    logging.debug("A total of {0} allowed ids loaded.".format(len(allowed)))

    newfrgfile = frgfile.replace(".frg", ".filtered.frg")
    fp = open(frgfile)
    fw = open(newfrgfile, "w")

    nfrags, discarded_frags = 0, 0
    nmates, discarded_mates = 0, 0
    for rec in iter_records(fp):
        if rec.type == "FRG":
            readname = rec.get_field("acc")
            readname = readname.rstrip("ab")
            nfrags += 1
            if readname not in allowed:
                discarded_frags += 1
                continue
        if rec.type == "LKG":
            readname = rec.get_field("frg")
            readname = readname.rstrip("ab")
            nmates += 1
            if readname not in allowed:
                discarded_mates += 1
                continue
        print(rec, file=fw)

    # Print out a summary
    survived_frags = nfrags - discarded_frags
    survived_mates = nmates - discarded_mates
    print(
        "Survived fragments: {0}".format(percentage(survived_frags, nfrags)),
        file=sys.stderr,
    )
    print(
        "Survived mates: {0}".format(percentage(survived_mates, nmates)),
        file=sys.stderr,
    )


def frg(args):
    """
    %prog frg frgfile

    Extract FASTA sequences from frg reads.
    """
    p = OptionParser(frg.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (frgfile,) = args
    fastafile = frgfile.rsplit(".", 1)[0] + ".fasta"
    fp = open(frgfile)
    fw = open(fastafile, "w")

    for rec in iter_records(fp):
        if rec.type != "FRG":
            continue
        id = rec.get_field("acc")
        seq = rec.get_field("seq")
        s = SeqRecord(Seq(seq), id=id, description="")
        SeqIO.write([s], fw, "fasta")

    fw.close()


def asm(args):
    """
    %prog asm asmfile

    Extract FASTA sequences from asm reads.
    """
    p = OptionParser(asm.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (asmfile,) = args
    prefix = asmfile.rsplit(".", 1)[0]
    ctgfastafile = prefix + ".ctg.fasta"
    scffastafile = prefix + ".scf.fasta"
    fp = open(asmfile)
    ctgfw = open(ctgfastafile, "w")
    scffw = open(scffastafile, "w")

    for rec in iter_records(fp):
        type = rec.type
        if type == "CCO":
            fw = ctgfw
            pp = "ctg"
        elif type == "SCF":
            fw = scffw
            pp = "scf"
        else:
            continue

        id = rec.get_field("acc")
        id = id.translate(None, "()").split(",")[0]
        seq = rec.get_field("cns").translate(None, "-")
        s = SeqRecord(Seq(seq), id=pp + id, description="")
        SeqIO.write([s], fw, "fasta")
        fw.flush()

    fw.close()


def count(args):
    """
    %prog count frgfile

    Count each type of messages
    """
    p = OptionParser(count.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (frgfile,) = args
    fp = open(frgfile)

    counts = defaultdict(int)
    for rec in iter_records(fp):
        counts[rec.type] += 1

    for type, cnt in sorted(counts.items()):
        print("{0}: {1}".format(type, cnt), file=sys.stderr)


if __name__ == "__main__":
    main()
