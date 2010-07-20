import sys, random, subprocess, os, re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIStandalone

#transform a fasta file generated from MetaSim. Change all BW sequences to FW sequences.
oldFile = sys.argv[1]
newFile = sys.argv[2]

handle = open(oldFile)
seqs = [seq_record for seq_record in SeqIO.parse(handle, 'fasta')]
handle.close()
num = len(seqs)

newSeqs = []
for i in range(num):
    if re.search(',bw,0-\w*;KEY=', seqs[i].description) or re.search(';KEY=[\w\.]*,fw,0-', seqs[i].description):
        continue;
    if re.search('bw', seqs[i].description):
        newSeqs.append(SeqRecord(seq = seqs[i].seq.reverse_complement(), \
                     id = seqs[i].id, \
                     description = seqs[i].description))
    else:
        newSeqs.append(seqs[i])

handle = open(newFile, 'w')
SeqIO.write(newSeqs, handle, 'fasta')
handle.close()