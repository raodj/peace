import os, sys

# Program is called with the following parameters
# argv[1] = geneFile
# argv[2] = numOf454TypeESTs
# argv[3] = "454" (literal)
# argv[5] = msOutputFile
# argv[6] = estFile454

# Program options
geneFile = sys.argv[1]
metasimType = sys.argv[3]

if len(sys.argv) > 4:
	mstFile = sys.argv[4]
	estFile = sys.argv[5]
else:
	mstFile = '.'
	estFile = '.'

# Run metasim with geneFile and numOf454TypeESTs to generate
os.system('python MetasimSetup%s.py %s %s' %(metasimType, geneFile, sys.argv[2]))
estFN = "%s-%s" %(geneFile.split('.')[0], metasimType)

#transform bw to fw
estFN_new = estFN + "_new"
os.system('python 454BwtoFw.py %s.fa %s.fa' %(estFN, estFN_new))
os.system('mv %s.fa %s.fa' %(estFN_new, estFN))

# Run the pipeline
os.system('python pipeline.py -e %s.fa -p -c 1' %(estFN)) 

# Clean up output
os.system('mv peace_estsim_'+estFN+'/peace_estsim_'+estFN+'.mst '+mstFile)
os.system('mv peace_estsim_'+estFN+'/peace_estsim_'+estFN+'.out clusters.out')
os.system('mv peace_estsim_'+estFN+'/estsim_'+estFN+'_fmt.fa '+estFile)
os.system('rm -rf peace_estsim_'+estFN+'/')

#os.system('rm peace_scr*.job')

# MST file will be in file mstFile
# ESTs will be in file estFile
# Clusters will be in file clusters.out
