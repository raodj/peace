#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re
import sys
import random

blocks = [""]
reverse_blocks = [""]
def getBlock(i, length):
    """Get block[abs(i)] fo blocks (generating new content if i > len(blocks),
    and reverse complement it if i < 0"""
    assert(i != 0)
    if (abs(i) >= len(blocks)):
        new_blocks = ["".join([random.choice(('A', 'C', 'G', 'T')) for x in range(length)]) for j in range(abs(i) + 1 - len(blocks))]
        new_reverse = [str(Seq(b, IUPAC.unambiguous_dna).reverse_complement()) for b in new_blocks]
        blocks.extend(new_blocks)
        reverse_blocks.extend(new_reverse)
    return blocks[i] if i > 0 else reverse_blocks[-i]

def main(argv):
    length = int(argv[1])
    in_file = argv[2]
    out_file = argv[3]

    ifp = open(in_file)
    ofp = open(out_file, "w")

    r = re.compile(",\s*")

    results = [(Seq("".join([getBlock(i,length) for i in [int(j) for j in r.split(line.rstrip())]]), IUPAC.unambiguous_dna), ",".join(r.split(line.rstrip()))) for line in ifp if line != '\n' and line[0] != '#']

    
    id = 0
    records = []
    for s,line in results:
        records.append(SeqRecord(s, id = "s" + str(id), description = line))
        id = id+1

    SeqIO.write(records, ofp, "fasta");

main(sys.argv)

