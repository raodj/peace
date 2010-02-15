#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re
import sys
import random

def randomSeq(length):
    """Return a random sequence of length l, bases choosen
    independetly from a uniform distribution"""
    return Seq("".join([random.choice(('A', 'C', 'G', 'T')) for x in range(length)]), IUPAC.unambiguous_dna);

def generateBlocks(n, length):
    """Generate n blocks of size length"""
    return [randomSeq(length) for x in range(0, length)]

def main(argv):
    in_file = argv[1]
    out_file = argv[2]

    ifp = open(in_file)
    ofp = open(out_file, "w")

    length = int(ifp.readline())
    n = int(ifp.readline())
    blocks = [""] + generateBlocks(n, length)

    r = re.compile(",\s*")
    results = [(Seq("".join([str(blocks[i]) if i > 0 else str(blocks[-i].reverse_complement()) for i in [int(j) for j in r.split(line.rstrip())]]), IUPAC.unambiguous_dna), ",".join(r.split(line.rstrip()))) for line in ifp if line != '\n' and line[0] != '#']

    
    id = 0
    records = []
    for (s,line) in results:
        records.append(SeqRecord(s, id = "s" + str(id), description = line))
        id = id+1

    SeqIO.write(records, ofp, "fasta");

main(sys.argv)

