#!/usr/bin/python

#--------------------------------------------------------------------
#
# This file is part of PEACE.
# 
# PEACE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PEACE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
# 
# Miami University makes no representations or warranties about the
# suitability of the software, either express or implied, including
# but not limited to the implied warranties of merchantability,
# fitness for a particular purpose, or non-infringement.  Miami
# University shall not be liable for any damages suffered by licensee
# as a result of using, result of using, modifying or distributing
# this software or its derivatives.
#
# By using or copying this Software, Licensee agrees to abide by the
# intellectual property laws, and all other applicable laws of the
# U.S., and the terms of GNU General Public License (version 3).
#
# Authors:   John Karro                    karroje@muohio.edu
#            Dhananjai M. Rao              raodm@muohio.edu
#
#---------------------------------------------------------------------

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
    if (len(argv) != 4):
        sys.stderr.write("Regression test generator.\n");
        sys.stderr.write("Copyright (C) Miami University, 2010-\n\n");
        sys.stderr.write("Number of parameters mismatched.\n");
        sys.stderr.write("Usage: <Block size (in nt)> <input template file name> <output fasta file name>\n");
        exit(1);
        
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

