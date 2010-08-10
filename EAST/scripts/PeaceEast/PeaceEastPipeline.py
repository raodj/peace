import sys, random, subprocess, os, re, time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def createJobFile(scriptName, proc, invocs, loadModules=False):
    if proc % 2 == 0:
        ppn = 2
        nodes = proc / 2
    else:
        ppn = 1
        nodes = proc

    if nodes == 1:
        ppn = 1
        nodes = 1

    # create the job file
    sFile = open(scriptName+'.job', 'w')

    sFile.write('#!/bin/bash -l\n')
    sFile.write('#PBS -N '+scriptName+'\n')

    sFile.write('#PBS -l nodes=%d:ppn=%d\n' %(nodes, ppn))
    sFile.write('#PBS -l walltime=10:00:00\n')
    sFile.write('#PBS -j oe\n')
    if nodes == 1 and ppn == 1:
        sFile.write('#PBS -q serial_restr\n')


    if loadModules:
        sFile.write('module load java\n')
    sFile.write('export VIADEV_DEFAULT_RETRY_COUNT=7\n')
    sFile.write('export VIADEV_DEFAULT_TIME_OUT=21\n')

    # change to the correct directory
    sFile.write('cd '+os.getcwd()+'\n')

    for invoc in invocs:
        sFile.write(invoc+'\n')
        
    sFile.close()

    return scriptName+'.job'
# end createJobFile

def runPeace(estFile):
    peaceOutputFile = 'peace_'+estFile
    peaceInvoc = 'time mpiexec ./peace '
    peaceInvoc+='--estFile '+estFile+' --output '+peaceOutputFile+'.out --output-mst-file '+peaceOutputFile+'.mst'
    jobFile = createJobFile('peace_scr_%s' %(estFile), 4, [peaceInvoc])
    peaceJobID = subprocess.Popen(['qsub', jobFile], stdout=subprocess.PIPE).communicate()[0].split(".")[0]
    print peaceJobID
	# Wait for script jobs to finish by checking for the output file
    while not os.path.exists('peace_scr_%s.o%s' %(estFile, peaceJobID)):
        time.sleep(10)
    print "Peace job finished"
    return peaceOutputFile+'.mst'

def runPeaceEast(estFile, options):
    peaceOutputFile = 'peace_'+estFile
    mstFile = peaceOutputFile+'.mst'
    eastConsensusFile = '%s.contigs' %(estFile)
    eastSingletonFile = '%s.singlets' %(estFile)
    eastNum = '%s.num' %(estFile)
    eastInfo = '%s.info' %(estFile)
    
    eastInvoc = []
    eastInvoc.append('mpiexec ./peace --estFile '+estFile+' --output '+peaceOutputFile+'.out --output-mst-file '+peaceOutputFile+'.mst')
    eastInvoc.append('./east %s %s %s %s %s %s > %s' %(options,estFile,mstFile,eastConsensusFile,eastSingletonFile,eastNum,eastInfo))
    
    jobFile = createJobFile('east_scr_%s' %(estFile), 1, eastInvoc)
    eastJobID = subprocess.Popen(['qsub', jobFile], stdout=subprocess.PIPE).communicate()[0].split(".")[0]
    print eastJobID
    
    return eastJobID

def runParallel(aceOutFlag, options, estFileName):
    #do all the opeation in the temporary directory tmpDir
    if os.path.exists('tmpDir'):
        subprocess.Popen(["rm", "-rf", "tmpDir"]).wait()
    os.mkdir("tmpDir")
    subprocess.Popen(['cp', estFileName, 'tmpDir/']).wait()
    subprocess.Popen(['cp', 'peace', 'tmpDir/']).wait()
    subprocess.Popen(['cp', 'SepReverseCluster', 'tmpDir/']).wait()
    subprocess.Popen(['cp', 'east', 'tmpDir/']).wait()
    os.chdir('tmpDir')
    
    #start main logic
    mstFile = runPeace(estFileName)
    
    #call the program to separate clusters(and make reverse complements) into different files by putting all the ESTs in one cluster into one file
    numCluster = int(subprocess.Popen(['./SepReverseCluster', estFileName, mstFile], stdout=subprocess.PIPE).communicate()[0].split("\n")[0])

    # submit all the subjobs, each job for each cluster
    eastJobIDs=[]
    for i in range(numCluster):
        fileName = estFileName + '.' + str(i)
        jobId = runPeaceEast(fileName, options)
        eastJobIDs.append(jobId)
        
    # wait  all the jobs finished
    # Wait for script jobs to finish by checking for the output file
    for i in range(numCluster):
        fileName = estFileName + '.' + str(i)
        while not (os.path.exists('east_scr_%s.o' %(fileName) +eastJobIDs[i]) or os.path.exists('east_scr_%s.e' %(fileName) +eastJobIDs[i])):
            time.sleep(10)
    
    #write the output contig file and singleton file
    outContigFile = estFileName + ".east.contigs"
    outSingletFile = estFileName + ".east.singlets"
    outAceFile = estFileName + ".east.ace"
    outTmpAceFile = estFileName + ".east.ace.tmp"
    handleC = open(outContigFile, 'w')
    handleS = open(outSingletFile, 'w')
    if aceOutFlag==1:
        handleATmp = open(outTmpAceFile, 'w')
    idxC = 1
    aceNumContig = 0
    aceNumReads = 0
    for i in range(numCluster):
        curContigFile = estFileName + '.' + str(i) + '.contigs'
        curSingFile = estFileName + '.' + str(i) + '.singlets'
        curAceFile = estFileName + '.' + str(i) + '.contigs.ace'
        
        handle = open(curContigFile)
        contigs = [seq_record for seq_record in SeqIO.parse(handle, 'fasta')]
        handle.close()
        for i in range(len(contigs)): #write these contigs into the output contig file
            handleC.write('>Contig'+str(idxC)+'\n')
            handleC.write(str(contigs[i].seq)+'\n')
            idxC += 1

        
        handle = open(curSingFile)
        singletons = [seq_record for seq_record in SeqIO.parse(handle, 'fasta')]
        SeqIO.write(singletons, handleS, 'fasta') #write these contigs into the output singleton file
        handle.close()
        
        if aceOutFlag==1:
            handle = open(curAceFile)
            firstLine = handle.readline().strip().split(' ') #first line: AF 1 233
            aceNumContig += int(firstLine[1])
            aceNumReads += int(firstLine[2])
            for line in handle:
                handleATmp.write(line)
            handle.close()
        
    handleC.close()
    handleS.close()
    if aceOutFlag==1:
        handleATmp.close()
        #write to outAceFile
        handleATmp = open(outTmpAceFile)
        handleA = open(outAceFile, 'w')
        handleA.write('AF ' + str(aceNumContig) + ' ' + str(aceNumReads) + '\n')
        for line in handleATmp:
            handleA.write(line)
        handleATmp.close()
        handleA.close()
        subprocess.Popen(['mv', outAceFile, '../']).wait()
    subprocess.Popen(['mv', outContigFile, '../']).wait()
    subprocess.Popen(['mv', outSingletFile, '../']).wait()

def runSerial(aceOutFlag, options, estFileName):
    if os.path.exists('tmpDir'):
        subprocess.Popen(["rm", "-rf", "tmpDir"]).wait()
    os.mkdir("tmpDir")
    subprocess.Popen(['cp', estFileName, 'tmpDir/']).wait()
    subprocess.Popen(['cp', 'peace', 'tmpDir/']).wait()
    subprocess.Popen(['cp', 'ProcessRC', 'tmpDir/']).wait()
    subprocess.Popen(['cp', 'east', 'tmpDir/']).wait()
    os.chdir('tmpDir')

    peaceOutputFile = 'peace_'+estFileName
    mstFile = peaceOutputFile+'.mst'
    eastConsensusFile = '%s.contigs' %(estFileName)
    eastSingletonFile = '%s.singlets' %(estFileName)
    eastNum = '%s.num' %(estFileName)
    eastInfo = '%s.info' %(estFileName)
    
    peaceCommand = 'mpiexec ./peace --estFile '+estFileName+' --output '+peaceOutputFile+'.out --output-mst-file '+peaceOutputFile+'.mst'
    subprocess.Popen(peaceCommand, shell=True).wait()
    subprocess.Popen(['./ProcessRC', estFileName, mstFile, '%s_tmp' %(estFileName)]).wait()
    subprocess.Popen(['mv', '%s_tmp' %(estFileName), estFileName]).wait()
    eastCommand = './east %s %s %s %s %s %s > %s' %(options,estFileName,mstFile,eastConsensusFile,eastSingletonFile,eastNum,eastInfo)
    print eastCommand
    subprocess.Popen(eastCommand, shell=True).wait()

    outContigFile = estFileName + ".east.contigs"
    outSingletFile = estFileName + ".east.singlets"
    outAceFile = estFileName + ".east.ace"

    if aceOutFlag==1:
        subprocess.Popen('mv %s.contigs.ace ../%s' %(estFileName, outAceFile), shell=True).wait()
    subprocess.Popen('mv %s.contigs ../%s' %(estFileName, outContigFile), shell=True).wait()
    subprocess.Popen('mv %s.singlets ../%s' %(estFileName, outSingletFile), shell=True).wait()


def main():
    #command to run the pipeline: python PeaceEastPipeline.py <estFile> -s|-p options
    estFileName = sys.argv[1]
    argc = len(sys.argv)
    if argc<3:
        print "wrong options"
    else: 
        aceOut = 0 #not output ace file
        marker = 3
        options = ""
        while marker<argc:
            switch = sys.argv[marker]
            marker += 1
            if (switch=="-ace") or (switch=="-OUTPUT_ACE"):
                aceOut = 1
            options += " " + switch
            
        if sys.argv[2]=="-s":
            runSerial(aceOut, options, estFileName)
        if sys.argv[2]=="-p":
            runParallel(aceOut, options, estFileName)

if __name__ == "__main__":
	main()  