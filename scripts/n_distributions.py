# Requires biopython to run. (module load biopython)
# Usage:
# python n_distributions.py [FastA filename] [--normalize]
# --normalize is an optional flag that will instruct
# the script to treat any non-"ACGT" character as an N
# (including lowercase letters and Xs).  This emulates the
# way that Peace handles such characters.  Default is false.

import sys, re
from Bio import SeqIO


normalize = False
if len(sys.argv) > 2 and sys.argv[2].strip() == "--normalize":
	normalize = True
numBases = 0
numN = 0
lengthDist = dict()
handle = open(sys.argv[1])
for seq_record in SeqIO.parse(handle, "fasta"):
	seq = seq_record.seq
	seqStr = str(seq)
	if normalize:
		# Normalize sequence
		seqStr = re.sub('[^ACGT]', 'N', seqStr)
	# Get the basic counts
	numBases += len(seqStr)
	numN += seqStr.count('N')
	# Add to the length distribution
	results = re.findall('N{1,}', seqStr)
	for pattern_match in results:
		matchLen = len(pattern_match)
		try:
			lengthDist[matchLen] = lengthDist[matchLen]+1
		except KeyError:
			lengthDist[matchLen] = 1
	
handle.close()
nPercent = float(numN) / float(numBases)
print nPercent
print lengthDist
