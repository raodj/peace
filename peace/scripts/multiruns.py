import os, sys

outputFile = open('multiRunOutput.txt', 'w')
outputFile.write('y_wcd\t\ty_peace\tr_wcd\t\tr_peace\n')
outputFile.close()

for i in range(0, int(sys.argv[1])):
        # Run the pipeline
	os.system('python pipeline.py all_zf_cdnas.fa zf.dig 1 sbe 0.010000 10 0.040000 1 1 10 0 0 0 -p --analyzer d2 --clusterMaker mst --estIdx 0 --frame 100 --word 6 -c 60 -w -c')
	# Do processing of output
	f = open('wcd_peace_estsim_all_zf_cdnas.fa_zf.dig[1_sbe_0.010000_10_0.040000_1_1_10_0_0_0_][--analyzer_d2_--clusterMaker_mst_--estIdx_0_--frame_100_--word_6_]/analysis.txt', 'r')
	younden = []
	rand = []
	for line in f:
		if line.split(':')[0] == "Younden Index":
			print line
			younden.append(line.split(':')[1].strip())
		elif line.split(':')[0] == "Uncorrected Rand Index":
			print line
			rand.append(line.split(':')[1].strip())
	f.close()
	outputFile = open('multiRunOutput.txt', 'a')
	outputFile.write('%s\t%s\t%s\t%s\n' %(younden[1], younden[0], rand[1], rand[0]))
	# Clean up output
	os.system('rm -rf wcd_peace_estsim_all_zf_cdnas.fa_zf.dig[1_sbe_0.010000_10_0.040000_1_1_10_0_0_0_][--analyzer_d2_--clusterMaker_mst_--estIdx_0_--frame_100_--word_6_]/')
	outputFile.close()
