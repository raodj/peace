import sys, subprocess, os, re, time, random
from Bio import SeqIO
from Bio.Blast import NCBIStandalone

def main(coverage, uniqueId, num_genes):
    geneSource = "all_zf_cdnas.reduced.fa"
    minGeneLen = 5000
    maxGeneLen = 8000
    errorRate = str(0.04)
    #coverage = int(sys.argv[1])
    #uniqueId = sys.argv[2]
    geneFile = "multiGeneFile_" + uniqueId
    meanLen = choose_genes(geneSource, minGeneLen, maxGeneLen, num_genes, geneFile)
    numOfEst = str(int(coverage*meanLen/800))
    estFile = "estFile_" + uniqueId
    mstFile = "mstFile_" + uniqueId #we won't use this file. We'll call peace to generate it because we need collect memory information for PEACE+EAST
    
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
    genEST(peaceGeneFile, numOfEst, errorRate, mstFile, estFile)
    
    eastId = runPeaceEast(estFile, mstFile)
    cap3Id = runCap3(estFile)
    tgiclId = runTgicl(estFile)
    velvetId = runVelvet(estFile)
    miraId = runMira(estFile)
    
    print 'wait on EAST Id: %s' %(eastId)
    wait_on_job(eastId)
    print 'wait on CAP3 Id: %s' %(cap3Id)
    wait_on_job(cap3Id)
    print 'wait on TGICL Id: %s' %(tgiclId)
    wait_on_job(tgiclId)
    print 'wait on VELVET Id: %s' %(velvetId)
    wait_on_job(velvetId)
    print 'wait on MIRA Id: %s' %(miraId)
    wait_on_job(miraId)
    writeInfo(eastId, cap3Id, tgiclId, velvetId, miraId)   
    
def wait_on_job(jobId):
    command = 'qstat %s >/dev/null 2> /dev/null' %(jobId)
    while os.system(command)==0: #if it's running, return 0; else return a positive interger (exit status, 39168)
        time.sleep(10)

def writeInfo(eastId, cap3Id, tgiclId, velvetId, miraId):
    outHandle = open('memory.out', 'w') 
    outHandle.write("EAST(kb)\tCAP3(kb)\tTGICL(kb)\tVELVET(kb)\tMIRA(kb)\n")
    outHandle.write('%s\t%s\t%s\t%s\t%s\n' %(getMemoryInfo(eastId), getMemoryInfo(cap3Id), getMemoryInfo(tgiclId), getMemoryInfo(velvetId), getMemoryInfo(miraId)))
    outHandle.close()

def getMemoryInfo(jobId):
    fileName = '/usr/local/torque/current/var/spool/torque/server_priv/accounting/%s' %(time.strftime("%Y%m%d", time.gmtime()))
    mem = subprocess.Popen(['grep', jobId, fileName], stdout=subprocess.PIPE).communicate()[0].strip()
    r = re.match('[\d\D]*resources_used.mem=((\d)+)kb[\d\D]*', mem)
    return r.group(1)

def choose_genes(geneSource, minLen, maxLen, num, ofile):
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

def genEST(geneFile, numOfEst, errorRate, mstFile, estFile):
#call Peace pipeline to generate est and MST files.
    baseDir = os.getcwd()
    os.chdir(baseDir+"/peace")
    #python wrapper.py $geneFile $numOfEst $errorRate $mstFile $estFile
    subprocess.Popen(['python', 'wrapper.py', geneFile, numOfEst, errorRate, mstFile, estFile]).wait()
    subprocess.Popen(['mv', estFile, '..']).wait()
    os.chdir(baseDir)
   
def createJobFile(jobFile, cmd):
    outHandle = open(jobFile, 'w') 
    outHandle.write("#!/bin/bash -l\n")
    outHandle.write("#PBS -N MemoryFootprint\n")
    outHandle.write("#PBS -l nodes=1:ppn=1 \n")
    outHandle.write("#PBS -l walltime=10:0:0 \n")
    outHandle.write("#PBS -m abe\n")
    outHandle.write("#PBS -V\n")
    #outHandle.write("cd $PBS_O_WORKDIR")
    outHandle.write(cmd)
    outHandle.close()

def runPeaceEast(estFile, mstFile):
    subprocess.Popen(['cp', estFile, 'east/'])
    os.chdir('east/')
    cmd = 'cd %s\n' %(os.getcwd())
    cmd = cmd + 'mpiexec ./peace --estFile %s --output-mst-file %s\n' %(estFile, mstFile)
    cmd = cmd + './Main %s %s %s.eastConsensusFile %s.eastSingletonFile %s.eastNum\n' %(estFile, mstFile, estFile, estFile, estFile)
    jobFileName = '%s.east.job' %(estFile)
    createJobFile(jobFileName, cmd)
    l = re.split('\.', subprocess.Popen(['qsub', jobFileName], stdout=subprocess.PIPE).communicate()[0].strip())
    os.chdir('..')
    return l[0]
    
def runCap3(estFile):
    subprocess.Popen(['cp', estFile, 'cap3/'])
    os.chdir('cap3/')
    cmd = 'cd %s\n' %(os.getcwd())
    cmd = cmd + 'cap3 %s' %(estFile)
    jobFileName = '%s.cap.job' %(estFile)
    createJobFile(jobFileName, cmd)
    l = re.split('\.', subprocess.Popen(['qsub', jobFileName], stdout=subprocess.PIPE).communicate()[0].strip())
    os.chdir('..')
    return l[0]
    
def runTgicl(estFile):
    subprocess.Popen(['cp', estFile, 'tgicl/'])
    os.chdir('tgicl/')
    cmd = 'cd %s\n' %(os.getcwd())
    cmd = cmd + 'tgicl %s' %(estFile)
    jobFileName = '%s.tgicl.job' %(estFile)
    createJobFile(jobFileName, cmd)
    l = re.split('\.', subprocess.Popen(['qsub', jobFileName], stdout=subprocess.PIPE).communicate()[0].strip())
    os.chdir('..')
    return l[0]

def runVelvet(estFile):    
    subprocess.Popen(['cp', estFile, 'velvet/'])
    os.chdir('velvet/')
    cmd = 'cd %s\n' %(os.getcwd())
    cmd = cmd + 'velveth . 21 -fasta -long %s\n' %(estFile)
    cmd = cmd + "velvetg . -cov_cutoff auto\n"
    jobFileName = '%s.velvet.job' %(estFile)
    createJobFile(jobFileName, cmd)
    l = re.split('\.', subprocess.Popen(['qsub', jobFileName], stdout=subprocess.PIPE).communicate()[0].strip())
    os.chdir('..')
    return l[0]

def runMira(estFile):
    subprocess.Popen(['cp', estFile, 'mira/msd_in.sanger.fasta'])
    os.chdir('mira/')
    cmd = 'cd %s\n' %(os.getcwd())
    cmd = cmd + 'mira -fasta -project=msd -job=denovo,est,sanger SANGER_SETTINGS -LR:wqf=no:mxti=no -AS:epoq=no\n'
    jobFileName = '%s.mira.job' %(estFile)
    createJobFile(jobFileName, cmd)
    l = re.split('\.', subprocess.Popen(['qsub', jobFileName], stdout=subprocess.PIPE).communicate()[0].strip())
    os.chdir('..')
    return l[0]
    
if __name__ == "__main__":
	main(30, '1', 1)  