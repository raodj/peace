# Setup for Illumina

import os, sys

os.system('MetaSim cmd -r %s -c -m -g errormodel-62bp.mconf %s' %(sys.argv[2], sys.argv[1]))
f1 = open('%s-Empirical.fna' %(sys.argv[1].split('.')[0]))
f2 = open('%s-Illumina.fa' %(sys.argv[1].split('.')[0]), 'w')
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
os.system("rm %s-Empirical.fna" %(sys.argv[1].split('.')[0]))
f2.close()
