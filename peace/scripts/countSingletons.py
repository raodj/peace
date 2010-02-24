# Usage: python countSingletons.py [FastA filename]

import sys

f1 = open(sys.argv[1], 'r')
count = 0
estCount = 0
for line in f1:
	if line[0] == 'C':
		if estCount == 1:
			count += 1
		estCount = 0
	elif line != "\n":
		estCount += 1
if estCount == 1:
	count+=1

print count
