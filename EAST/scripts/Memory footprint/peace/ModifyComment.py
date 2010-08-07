import sys, random, subprocess, os, re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIStandalone

#transform a fasta file generated from MetaSim. Change all BW sequences to FW sequences.
oldFile = sys.argv[1]
newFile = oldFile + ".tmp"

handle = open(oldFile)
seqs = [seq_record for seq_record in SeqIO.parse(handle, 'fasta')]
handle.close()
num = len(seqs)

newSeqs = []
for i in range(num):
        newSeqs.append(SeqRecord(seq = seqs[i].seq, \
                     id = "%s_%d" %(seqs[i].id, i), \
                     description = "%s_%d" %(seqs[i].description, i)))

handle = open(newFile, 'w')
SeqIO.write(newSeqs, handle, 'fasta')
handle.close()

subprocess.Popen(['mv', newFile, oldFile]).wait()