# Script to annotate an EST file with the source genes
# as part of each EST's identifier

import sys

estf = open(sys.argv[1], 'r') # EST file
linkf = open(sys.argv[2], 'r') # file linking genes to EST names
newestf = open(sys.argv[3], 'w') # new EST file

# set up dictionary mapping EST names to proper identifiers
linkDict = dict()
for line in linkf:
	entries = line.split()
	ident = '>'+entries[0]+'_'+entries[1]+'_'+entries[2]
	linkDict[entries[4]] = ident
linkf.close()

# read estf and write to newestf
for line in estf:
	if line[0] == '>':
		# identifier line so find the identifier in dict
		newestf.write(linkDict[line[1:].strip()]+'\n')
	else:
		# EST data so write the data
		newestf.write(line)
estf.close()
newestf.close()

# done
