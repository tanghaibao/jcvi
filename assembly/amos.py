#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Based on AMOS `amos-3.0.0/src/PythonModules`

AMOS message reader/writer
Contributed by Paul Harrison
"""

import re
import sys

from optparse import OptionParser
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcvi.apps.base import ActionDispatcher, debug
debug()


class Message (list):
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
        self.fields = defaultdict(str) 
    
    def __str__(self):
        result = '{' + self.type + '\n'
        for key, value in sorted(self.fields.items()):
            result += key + ':'
            if '\n' in value:
                result += '\n' + value 
                if not value.endswith('\n'):
                    result += '\n'
                result += '.\n'
            else:
                result += value +'\n'        

        result += ''.join(str(sub_message) for sub_message in self) + '}\n'

        return result


_START = re.compile(r'^{([A-Z][A-Z][A-Z])\n$')
_MULTILINE_FIELD = re.compile(r'^([a-z][a-z][a-z]):\n$')
_FIELD = re.compile(r'^([a-z][a-z][a-z]):(.*)\n$')


def read_record(fp, first_line=None):
    """
    Read a record from a file of AMOS messages 
    
    On success returns a Message object     
    On end of file raises EOFError  
    On syntax error raises amos.Error
    """

    if first_line is None:
        first_line = fp.next()
    
    if not first_line:
        raise EOFError()
    
    match = _START.match(first_line)
    if not match:
        raise Exception('Bad start of message', first_line)
    
    type = match.group(1)
    message = Message(type)
    
    while True:
        row = fp.next()
        
        match = _MULTILINE_FIELD.match(row)
        if match:
            key = match.group(1)
            while True:
                row = fp.next()
                if row[0] == '.': break
                message.fields[key] += row.strip()
            continue

        match = _FIELD.match(row)
        if match:
            key, val = match.group(1), match.group(2)
            message.fields[key] = val 
            continue

        if row[0] == '}':
            break
    
        if row[0] == '{':
            message.append( read_record(fp, row) )
            continue
        
        raise Exception('Bad line', row)

    return message


def iter_records(file):
    """ Iterate over all the records in a file """

    while True:
        try:
            yield read_record(file)
        except EOFError:
            break


def main():

    actions = (
        ('frg', 'extract FASTA sequences from frg'),
        ('count', 'count each type of messages'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def frg(args):
    """
    %prog frg frgfile

    Example: pass a message stream from stdin to stdout
             print count of each message type to stderr
    """
    p = OptionParser(frg.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    frgfile, = args
    fastafile = frgfile.rsplit(".", 1)[0] + ".fasta"
    fp = open(frgfile)
    fw = open(fastafile, "w")

    for rec in iter_records(fp):
        if rec.type != "FRG": continue
        id = rec.fields["acc"]
        seq = rec.fields["seq"]
        s = SeqRecord(Seq(seq), id=id, description="")
        SeqIO.write([s], fw, "fasta")

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

    frgfile, = args
    fp = open(frgfile)

    counts = defaultdict(int)
    for rec in iter_records(fp):
        counts[rec.type] += 1 
    
    for type, cnt in sorted(counts.items()):
        print >> sys.stderr, '{0}: {1}'.format(type, cnt)


if __name__ == '__main__':  
    main()
