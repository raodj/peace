import sys, random, subprocess, os, re
from Bio import SeqIO
from Bio.Blast import NCBIStandalone

geneSource = sys.argv[1]
num_genes = int(sys.argv[2])
minGeneLen = int(sys.argv[3])
maxGeneLen = int(sys.argv[4])
errorRate = float(sys.argv[5])/100
coverage = int(sys.argv[6])
useQualityFile = sys.argv[7] #1-use; 0-not use.
outputFile = sys.argv[8];
uniqueId = sys.argv[9];
dataType = sys.argv[10];

def choose_genes(minLen, maxLen, num, ofile):
    handle = open(geneSource)
    genes = [seq_record for seq_record in SeqIO.parse(handle, 'fasta') if len(seq_record.seq)>=minLen and len(seq_record.seq)<=maxLen]
    selectedGenes=[] 
    totalLen = 0
    if num < len(genes):
        for i in range(num):
            idx = random.randint(0, len(genes)-1)
            selectedGenes.append(genes[idx])
            totalLen += len(str(genes[idx].seq))
            del genes[idx]
    else:
        selectedGenes = genes
        for i in range(len(genes)):
            totalLen += len(str(genes[i].seq))
    handle.close()
    
    #write the selected genes to a file
    out = open(ofile, "w")
    SeqIO.write(selectedGenes, out, 'fasta')
    out.close()
    return totalLen/num
    

"""read all the contigs in the input file, put it into the hashTable as needed.
The hashTable key=gene name, value=the longest contig corresponding to the gene"""
def parseContigs(fileName, hashTable, contigToGene, dbFile):
    geneEValue = dict() #store the smallest E-value for the gene
    for x in hashTable.keys(): #initialize geneEValue 
        geneEValue[x] = 1000

    handle = open(fileName)
    contigs = [seq_record for seq_record in SeqIO.parse(handle, 'fasta')]
    handle.close()
    num = len(contigs)
    
    tmpContig = 'tmpContig.fa.%s' %(uniqueId)
    for i in range(num):
        contig = str(contigs[i].seq).strip()
        if len(contig) < 300:
            contigToGene.append('-')
            continue
        
        tmpIn = open(tmpContig, 'w')
        tmpIn.write('>contig\n')
        tmpIn.write(contig)
        tmpIn.close()
        f = 'blastOut.' + uniqueId
        subprocess.Popen(['blastall', '-i', tmpContig, '-d', dbFile, '-p', 'blastn', '-o', f]).wait()

        result_handle = open(f)
        blast_parser = NCBIStandalone.BlastParser()
        blast_record = blast_parser.parse(result_handle)
        ali = blast_record.alignments
        result_handle.close()    
        if len(ali)<0:
            contigToGene.append('-')
            continue
        key = re.split('>|\s+', ali[0].title)[1] + '_' + uniqueId
        contigToGene.append(key)

        if (ali[0].hsps[0].expect <= geneEValue[key]): #ali[0] is the closest match with the smallest E-value
            #For those having same E-value, choose the one with the longest contig.
            #If the length of the gene is 10000, one consensus is 10000-base long and same as the gene, another one is 2000-base long and is part of the gene, they will have the same E-value 0.0
            #similarly, if the length of the gene is 10000, one consensus is 20000-base long and includes the gene, the E-value will also be 0.0
            #a very small e-value means that they contain a significant subsequence that is the same, but cannot tell if they are same.
            curLen = len(hashTable[key])
            if len(contig)>curLen:
                hashTable[key] = contig
        subprocess.Popen(["rm", f]).wait()
    subprocess.Popen(["rm", tmpContig]).wait()
    
def getCoverageRate(congtigFileName, contigToGene, lenOfEachGene, geneSeqs, allKeys):
    geneCoverage = dict() #key-geneName, value- percentage of coverage
    geneToCon = dict() #key- geneName, value- a list including the indices of all the contigs which correspond to the gene
    for i in range(len(allKeys)):
        geneToCon[allKeys[i]] = []
    for i in range(len(contigToGene)):
        if contigToGene[i]=='-':
            continue
        geneToCon[contigToGene[i]].append(i)
    
    handle = open(congtigFileName)
    contigs = [seq_record for seq_record in SeqIO.parse(handle, 'fasta')]
    handle.close()
    
    for key in allKeys:
        idxOfContigs = geneToCon[key]
        if len(idxOfContigs)==0:
            geneCoverage[key] = 1
            continue
            
        coveredPos = dict() #key- index of position, value- 1 indicates the position is covered by some contig
        for i in idxOfContigs:
            contig = str(contigs[i].seq).strip()
            if len(contig) == 0:
                continue
            gene = geneSeqs[key]
            command = 'java -jar SWSmart.jar %s %s' %(contig, gene)
            pos = os.popen(command, 'r').readline().strip().split(' ') # "startPos  size"
            startPos = int(pos[0])
            length = int(pos[1])
            if startPos == -1:
                continue
            for x in range(length):
                coveredPos[startPos+x] = 1
        geneCoverage[key] = 1-float(coveredPos.values().count(1))/lenOfEachGene[key]        
    return geneCoverage
    
#------------------------------------------------------------------------------------------------
#choose genes
geneFile = "multiGeneFile_" + uniqueId
meanLen = choose_genes(minGeneLen, maxGeneLen, num_genes, geneFile)

numOfEst = ''
if dataType == "sanger":
    numOfEst = str(int(coverage*meanLen/800))
elif dataType == "454":
    numOfEst = str(int(coverage*meanLen/230)*num_genes)
elif dataType == "illumina":
    numOfEst = str(int(coverage*meanLen/62)*num_genes)


baseDir = os.getcwd()

#peace wrapper only accepts gene file with one gene on one line rather than in multiple lines. 
#here we transform geneFile into the format accepted by the wrapper
source = open(geneFile)
tGenes = [seq_record for seq_record in SeqIO.parse(source, 'fasta')]
source.close()
peaceGeneFile = geneFile + '_Peace'
outPeace = open(peaceGeneFile, 'w')
for i in range(len(tGenes)):
    outPeace.write('>'+tGenes[i].description+'\n')
    outPeace.write(str(tGenes[i].seq)+'\n')
outPeace.close()

subprocess.Popen(['mv', peaceGeneFile, 'peace/']).wait()

usedTime = subprocess.Popen(['./singleTestForMulti.pl', peaceGeneFile, numOfEst, str(errorRate), uniqueId, useQualityFile, baseDir, dataType], stdout=subprocess.PIPE).communicate()[0].strip()

#The following code analyzes the result
estFileName = 'estFile_%s.fa' %(uniqueId)
eastConsensusFile = 'estFile_%s.fa.east.contigs' %(uniqueId)
eastSingletonFile = "estFile_%s.fa.east.singlets" %(uniqueId)
cap3ConsensusFile = "estFile_%s.fa.cap.contigs" %(uniqueId)
cap3SingletonFile = "estFile_%s.fa.cap.singlets" %(uniqueId)
tgiclConsensusFile = 'estFile_%s.fa.tgicl.contigs' %(uniqueId)
tgiclSingletonFile = "estFile_%s.fa.tgicl.singlets" %(uniqueId)
velvetConsensusFile = 'estFile_%s.fa.velvet.contigs' %(uniqueId)
miraConsensusFile = 'estFile_%s.fa.mira.contigs' %(uniqueId)
subprocess.Popen(['formatdb', '-i', geneFile, '-p', 'F']).wait() #generate blast db

#read out every gene in geneFile, put it  into different files named as the gene's description
handle = open(geneFile)
genes = [seq_record for seq_record in SeqIO.parse(handle, 'fasta')]
handle.close()
tmpGeneFiles = 'tmpGeneFiles' + '_' + uniqueId
if os.path.exists(tmpGeneFiles):
    subprocess.Popen(['rm', '-rf', tmpGeneFiles]).wait()
subprocess.Popen(['mkdir', tmpGeneFiles]).wait()
os.chdir(tmpGeneFiles)

num = len(genes)
geneToContigEAST = dict() #put the contig with the smallest E-value and  corresponding to this gene for EAST
contigToGeneEAST = []
geneToContigCAP = dict() #put the contig with the smallest E-value and  corresponding to this gene for CAP3
contigToGeneCAP = []
geneToContigTGICL = dict() #put the contig with the smallest E-value and  corresponding to this gene for TGICL
contigToGeneTGICL = []
geneToContigVELVET = dict() #put the contig with the smallest E-value and  corresponding to this gene for VELVET
contigToGeneVELVET = []
geneToContigMIRA = dict() #put the contig with the smallest E-value and  corresponding to this gene for MIRA
contigToGeneMIRA = []
geneLength = dict() #store the length of each gene for analysis
geneSeq = dict() #store the sequence of each gene
for i in range(num):
    desc = re.split('\s+', genes[i].description)[0] #get gene name
    desc = desc + '_' + uniqueId
    
    geneToContigEAST[desc] = ""
    geneToContigCAP[desc] = ""
    geneToContigTGICL[desc] = ""
    geneToContigVELVET[desc] = ""
    geneToContigMIRA[desc] = ""
    geneLength[desc] = len(genes[i].seq)
    geneSeq[desc] = genes[i].seq
    #write the gene to a file
    out = open(desc, "w")
    out.write(str(genes[i].seq)) #write the gene into the file named desc
    out.close()
os.chdir('..')

parseContigs(eastConsensusFile, geneToContigEAST, contigToGeneEAST, geneFile) #parse contigs generated by EAST
parseContigs(cap3ConsensusFile, geneToContigCAP, contigToGeneCAP, geneFile) #parse contigs generated by CAP3
parseContigs(tgiclConsensusFile, geneToContigTGICL, contigToGeneTGICL, geneFile) #parse contigs generated by TGICL
parseContigs(velvetConsensusFile, geneToContigVELVET, contigToGeneVELVET, geneFile) #parse contigs generated by VELVET
parseContigs(miraConsensusFile, geneToContigMIRA, contigToGeneMIRA, geneFile) #parse contigs generated by MIRA
keys = geneToContigEAST.keys()

covEAST = getCoverageRate(eastConsensusFile, contigToGeneEAST, geneLength, geneSeq, keys)
covCAP = getCoverageRate(cap3ConsensusFile, contigToGeneCAP, geneLength, geneSeq, keys)
covTGICL = getCoverageRate(tgiclConsensusFile, contigToGeneTGICL, geneLength, geneSeq, keys)
covVELVET = getCoverageRate(velvetConsensusFile, contigToGeneVELVET, geneLength, geneSeq, keys)
covMIRA = getCoverageRate(miraConsensusFile, contigToGeneMIRA, geneLength, geneSeq, keys)

finalOut = open(outputFile, 'w')
finalOut.write("geneName\tgeneLen\tCap3Time\tPeaceTime\tEastTime\tTgiclTime\tVelvetTime\tMiraTime\tCap3LenOfContig\tEastLenOfContig\tCap3AScore\tEastAScore\tCap3NumContigs\tEastNumContigs\tCap3NumSing\tEastNumSing\tnumEsts\tTgiclNumContigs\tVelvetNumContigs\tMiraNumContigs\tTgiclAScore\tVelvetAScore\tTgiclNumSing\tNonCovEAST\tNonCovCAP\tNonCovTGICL\tNonCovVELVET\tNonCovMIRA\tMiraAScore\n")

handle = open(cap3ConsensusFile)
cap3NumContigs = len([seq_record for seq_record in SeqIO.parse(handle, 'fasta')])
handle.close()
handle = open(eastConsensusFile)
eastNumContigs = len([seq_record for seq_record in SeqIO.parse(handle, 'fasta')])
handle.close()
handle = open(tgiclConsensusFile)
tgiclNumContigs = len([seq_record for seq_record in SeqIO.parse(handle, 'fasta')])
handle.close()
handle = open(velvetConsensusFile)
velvetNumContigs = len([seq_record for seq_record in SeqIO.parse(handle, 'fasta')])
handle.close()
handle = open(miraConsensusFile)
miraNumContigs = len([seq_record for seq_record in SeqIO.parse(handle, 'fasta')])
handle.close()

handle = open(cap3SingletonFile)
cap3NumSing = len([seq_record for seq_record in SeqIO.parse(handle, 'fasta')])
handle.close()
handle = open(eastSingletonFile)
eastNumSing = len([seq_record for seq_record in SeqIO.parse(handle, 'fasta')])
handle.close()
handle = open(tgiclSingletonFile)
tgiclNumSing = len([seq_record for seq_record in SeqIO.parse(handle, 'fasta')])
handle.close()
handle = open(estFileName)
numEsts = len([seq_record for seq_record in SeqIO.parse(handle, 'fasta')])
handle.close()

for x in keys:
    os.chdir(tmpGeneFiles)

    f1 = "east.out." + uniqueId
    out = open(f1, "w") #'east.out' is the default file name for eSTAssembly.ResultAnalysis
    out.write(geneToContigEAST[x]) 
    out.close()
    f1 = tmpGeneFiles + '/' + f1
    f2 = "cap3.out." + uniqueId
    out = open(f2, "w") #'cap3.out' is the default file name for eSTAssembly.ResultAnalysis
    out.write(geneToContigCAP[x]) 
    out.close()
    f2 = tmpGeneFiles + '/' + f2
    f3 = tmpGeneFiles + '/' + x
    f4 = 'analysis.out.' + uniqueId
    f5 = "tgicl.out." + uniqueId
    out = open(f5, "w") 
    out.write(geneToContigTGICL[x]) 
    out.close()
    f5 = tmpGeneFiles + '/' + f5
    f6 = "velvet.out." + uniqueId
    out = open(f6, "w") 
    out.write(geneToContigVELVET[x]) 
    out.close()
    f6 = tmpGeneFiles + '/' + f6
    f7 = "mira.out." + uniqueId
    out = open(f7, "w") 
    out.write(geneToContigMIRA[x]) 
    out.close()
    f7 = tmpGeneFiles + '/' + f7

    os.chdir('..')
    subprocess.Popen(['java', 'eSTAssembly.ResultAnalysis', '0', '0', '0', '0', 'analysis.properties', f3, f1, 'eastSingleton.out', f2, 'cap3Singleton.out', f4, f5, f6, f7]).wait()
    anaFile = open(f4)
    anaFile.readline() #skip the first two-line comments
    anaFile.readline()
    analysis = re.split('\t', anaFile.readline().strip())
    finalOut.write('%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n' %(x,geneLength[x], usedTime,len(geneToContigCAP[x]), len(geneToContigEAST[x]), analysis[6],analysis[7],cap3NumContigs,eastNumContigs,cap3NumSing,eastNumSing,numEsts,tgiclNumContigs,velvetNumContigs,miraNumContigs,analysis[12],analysis[13],tgiclNumSing,covEAST[x],covCAP[x],covTGICL[x],covVELVET[x],covMIRA[x],analysis[14]))

    os.remove(f4)
finalOut.write('DONE')
finalOut.close()    


#remove those intermediate files
subprocess.Popen(['rm', geneFile])
subprocess.Popen(['rm', eastConsensusFile])
subprocess.Popen(['rm', cap3ConsensusFile])
subprocess.Popen(['rm', tgiclConsensusFile])
subprocess.Popen(['rm', velvetConsensusFile])
subprocess.Popen(['rm', miraConsensusFile])
subprocess.Popen(['rm', eastSingletonFile])
subprocess.Popen(['rm', cap3SingletonFile])
subprocess.Popen(['rm', tgiclSingletonFile])
subprocess.Popen(['rm', estFileName])
subprocess.Popen(['rm', '-rf', tmpGeneFiles])

#for x in keys:
#    subprocess.Popen(['rm', x])
#remove the intermediate file generated by blast
f1 = geneFile + '.nin'
subprocess.Popen(['rm', f1]).wait()
f1 = geneFile + '.nsq'
subprocess.Popen(['rm', f1]).wait()
f1 = geneFile + '.nhr'
subprocess.Popen(['rm', f1]).wait()



