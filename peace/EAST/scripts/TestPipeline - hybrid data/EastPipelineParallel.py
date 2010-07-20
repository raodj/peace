import sys, os, time, subprocess

def wait_on_job(jobId):
    command = 'qstat %s >/dev/null 2> /dev/null' %(jobId)
    while os.system(command)==0: #if it's running, return 0; else return a positive interger (exit status, 39168)
        time.sleep(10)

if sys.argv[1] == '-r':
    batchFile = sys.argv[2]  
    outputFile = sys.argv[3]
    subprocess.Popen(['perl', 'EastPipelineParallel.pl', sys.argv[1], batchFile]).wait()
else:
    batchFile = sys.argv[1]        
    outputFile = sys.argv[2]
    subprocess.Popen(['perl', 'EastPipelineParallel.pl', batchFile]).wait()
    

handle = open(batchFile, 'r')
resultsFiles = []
for line in handle:
    tmp = line.split(' ')
    jobId = tmp[2].strip()
    wait_on_job(jobId)
    oFile = tmp[1]+'.o'+tmp[2]
    subprocess.Popen(['rm', oFile], stderr=subprocess.PIPE)
    oFile = tmp[1]+'.e'+tmp[2]
    subprocess.Popen(['rm', oFile], stderr=subprocess.PIPE)
    resultsFiles.append(tmp[4].strip())
handle.close()

'''
x=[0,1,2,3,4,5,6]
y=[20,30,40,50,60]
outputFile = sys.argv[2]
resultsFiles = []

for i in x:
    for j in y:
        for k in range(30):
            tmpStr = 'results.d/ng15.d/e%d.d/c%d.d/t%d.d/output.ng15_e%d_c%d_t%d' %(i,j,k+1,i,j,k+1)
            resultsFiles.append(tmpStr)
'''

#collect results
outHandle = open(outputFile, 'w') 
outHandle.write("No.\tTestNo.\tgeneName\tgeneLen\tErrorRate\tCoverage\tCap3LenOfContig\tEastLenOfContig\tCap3AScore\tEastAScore\tTgiclAScore\tVelvetAScore\tMiraAScore\n")
outputFile2 = outputFile + '.runtime'
outHandle2 = open(outputFile2, 'w')   
outHandle2.write("TestNo.\tNumOfGenes\tErrorRate\tCoverage\tCap3Time\tPeaceTime\tEastTime\tTgiclTime\tVelvetTime\tMiraTime\tCap3NumContigs\tEastNumContigs\tCap3NumOfSingletons\tEastNumOfSingletons\tNumOfESTs\tTgiclNumContigs\tTgiclNumOfSingletons\tVelvetNumContigs\tMiraNumContigs\n")

testNo = 1
num = 1
for x in resultsFiles:
    if (not os.path.exists(x)):
        continue
    tHandle = open(x)
    lineList = tHandle.readlines()
    print x
    if lineList.pop() == 'DONE':
        lineList.remove(lineList[0]) #remove the comment line
        lineInfo = lineList[0].split('\t')
        geneInfo = lineInfo[0].split('_') #format: TC301706_ng2_e0_c10_t1.1645
        numOfGenes = geneInfo[1][2:]
        errorRate = geneInfo[2][1:]
        coverage = geneInfo[3][1:]
        cap3Time = lineInfo[2]
        peaceTime = lineInfo[3]
        eastTime = lineInfo[4]
        tgiclTime = lineInfo[5]
        velvetTime = lineInfo[6]
        MiraTime = lineInfo[7]
        outHandle2.write('%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(testNo, numOfGenes, errorRate, coverage, cap3Time, peaceTime, eastTime, tgiclTime, velvetTime, MiraTime, lineInfo[12], lineInfo[13], lineInfo[14], lineInfo[15], lineInfo[16], lineInfo[17], lineInfo[22], lineInfo[18], lineInfo[19]))
        
        for x in range(len(lineList)):
            lineInfo = lineList[x].split('\t')
            outHandle.write('%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(num, testNo, lineInfo[0].split('_')[0], lineInfo[1], errorRate, coverage, lineInfo[8], lineInfo[9], lineInfo[10], lineInfo[11], lineInfo[20], lineInfo[21], lineInfo[23]))
            num += 1
        testNo += 1
        
    tHandle.close()
outHandle.close()
outHandle2.close()

