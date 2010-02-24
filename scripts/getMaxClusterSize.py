# Usage: python getMaxClusterSize.py [FastA filename]

import sys

f1 = open(sys.argv[1], 'r')
maxEstCount = 0
estCount = 0
for line in f1:
	if line[0] == 'C':
		if estCount > maxEstCount:
			maxEstCount = estCount
		if estCount > 100:
			print estCount
		estCount = 0
	elif line != "\n":
		estCount += 1

if line[0] == 'C':
	if estCount > maxEstCount:
		maxEstCount = estCount
	if estCount > 100:
		print estCount


print maxEstCount
