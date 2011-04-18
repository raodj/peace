# Setup for 454

# The program is called with the following opetions
# argv[1] = SourceGeneFile (input)
# argv[2] = numOf454ESTs to generate (input)

import os, sys

os.system('MetaSim cmd -r %s -c -4 %s' %(sys.argv[2], sys.argv[1]))
f1 = open('%s-454.fna' %(sys.argv[1].split('.')[0]))
f2 = open('%s-454.fa' %(sys.argv[1].split('.')[0]), 'w')
gNum = 0
cur = ""
for line in f1:
	if line[0] == '>':
		gLabel = line.strip().split("SOURCE_1=")[1]
		gLabel = gLabel.split(" ")[0]
		if cur != gLabel:
			cur = gLabel
			gNum+=1
		cLabel = "g%06d" %(gNum)
		f2.write(">%s_%s" %(cLabel, line[1:]))
	else:
		f2.write(line)
f1.close()
os.system("rm %s-454.fna" %(sys.argv[1].split('.')[0]))
f2.close()
