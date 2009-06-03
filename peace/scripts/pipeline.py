import os, sys, time, popen2
from optparse import OptionParser



class Pipeline:
	def __init__(self, specifyESTs, faFile, digFile, peace, wcd):
		self.specifyESTs = specifyESTs
		self.runningPeace = peace
		self.runningWcd = wcd
		if specifyESTs:
			self.estFile = faFile
		else:
			self.faFile = faFile
			self.digFile = digFile
		self.wcdMatrixOut = False
		self.proc = 1
	
	def setESTParams(self, estSimParams):
		self.estSimParams = estSimParams
		
	def setProcessors(self, proc):
		self.proc = proc
		
	def run(self):
		# Setup things
	
		# new code for specifying ESTs as an input file
		if self.specifyESTs:
			# copy file, so we save the original
			os.system('cp '+self.estFile+' estsim_'+self.estFile)
			estOutputFile = 'estsim_'+self.estFile[:-3]
		else:
			# set names of intermediate files, and strings used to invoke the programs
			estOutputFile = 'estsim_'+self.faFile+'_'+self.digFile+'['
			for x in self.estSimParams:
				estOutputFile+=x+'_'
			estOutputFile+=']'
	
			estSimInvoc = 'java -jar estsim.jar '+self.faFile+' '+self.digFile+' -o '+estOutputFile+'.fa'
			for x in self.estSimParams:
				estSimInvoc += ' '+x

	       	formatInvoc = 'python format.py '+estOutputFile+'.fa'
	
		dirName = estOutputFile
	
		if self.runningPeace:
		        peaceOutputFile = 'peace_'+estOutputFile
		        dirName = peaceOutputFile
			peaceInvoc = 'time mpiexec ./peace '
			peaceInvoc+='--estFile '+estOutputFile+'_fmt.fa --output '+peaceOutputFile+'.out --output-mst-file mstFile.mst'
	
		if self.runningWcd:
			dirName = 'wcd_'+estOutputFile
			wcdInvoc = 'time mpiexec ./wcd -b'
			if self.proc > 1:
				wcdInvoc+=' -N %d' %(self.proc)
			wcdInvoc += ' -o wcd_'+estOutputFile+'.txt '+estOutputFile+'_fmt.fa'
	
		# directory in which to store output and intermediate files
		if self.runningPeace and self.runningWcd:
			dirName = 'wcd_'+peaceOutputFile
	
		if self.specifyESTs:
			 # run format script only, not estsim
		       	formatInvoc = 'python formatNoSelect.py '+estOutputFile+'.fa'
			estSimInvocs = [formatInvoc]
			
		elif self.wcdMatrixOut:
			# need sorted order to produce matrix
			# if we are just generating ESTs, sort them too
			formatInvoc = 'python sort.py '+estOutputFile+'.fa'
			estSimInvocs = [estSimInvoc, formatInvoc]

		else:
			# normal case
			estSimInvocs = [estSimInvoc, formatInvoc]
			
		# create and run job(s)
		if self.runningPeace and self.runningWcd:
		        # need to create and run 3 jobs
		        # estsim will run first
		        self.createJobFile('estsim_script', 1, estSimInvocs, True)
		        p = popen2.popen2('qsub estsim_script.job')
			jobID = p[0].read().split('.')[0]
			# Wait for script jobs to finish by checking for the output file
			while not os.path.exists('estsim_script.o'+jobID):
		               	time.sleep(10)
		        # now create both peace and wcd job files, and run them in parallel
		        self.createJobFile('wcd_script', self.proc, [wcdInvoc])
		        self.createJobFile('peace_script', self.proc, [peaceInvoc])
		        #startTime = time.time()
			p = popen2.popen2('qsub peace_script.job')
			peaceJobID = p[0].read().split('.')[0]
			p = popen2.popen2('qsub wcd_script.job')
			wcdJobID = p[0].read().split('.')[0]
			print "JobID for peace_script: "+peaceJobID
			print "JobID for wcd_script: "+wcdJobID
			while not os.path.exists('peace_script.o'+peaceJobID) or not os.path.exists('wcd_script.o'+wcdJobID):
		                time.sleep(10)
		
		else:
			scriptName = 'estsim_script'
		        if self.runningPeace:
		                estSimInvocs.append(peaceInvoc)
		                scriptName = 'peace_script'
		        elif self.runningWcd:
		                estSimInvocs.append(wcdInvoc)
		                scriptName = 'wcd_script'
		        self.createJobFile(scriptName, self.proc, estSimInvocs, True)
		        startTime = time.time()
		        p = popen2.popen2('qsub '+scriptName+'.job')
		        jobID = p[0].read().split('.')[0]
		        print 'JobID for '+scriptName+': '+jobID
		        # Wait for script jobs to finish by checking for the output file
		        while not os.path.exists(scriptName+'.o'+jobID):
		                time.sleep(10)
		
		print "Jobs finished"
		
		# Organize output in a directory
		cwd = os.getcwd()
		entries = os.listdir(cwd)
		try:
			os.mkdir(dirName)
		except OSError:
			# Directory exists, make one with a unique identifier
			dirName = dirName+jobID
			os.mkdir(dirName)

		if not self.runningPeace and not self.runningWcd:
			os.system('mv '+estOutputFile+'.fa '+dirName+'/')
		
		else:
			os.system('mv '+estOutputFile+'_fmt.fa '+dirName+'/')
			os.system('mv '+estOutputFile+'_dict.txt '+dirName+'/')
		
		if self.runningPeace and self.runningWcd:
			os.system('mv estsim_script.o'+jobID+' '+dirName+'/')
			os.system('mv peace_script.o'+peaceJobID+' '+dirName+'/')
			os.system('mv wcd_script.o'+wcdJobID+' '+dirName+'/')
			os.system('rm peace_script.job')
			os.system('rm estsim_script.job')
			os.system('rm wcd_script.job')
		else:
			#os.system('rm '+scriptName+'.job')
			os.system('mv '+scriptName+'.o'+jobID+' '+dirName+'/')
			peaceJobID = jobID
			wcdJobID = jobID
		
		if self.runningPeace:
			os.system('mv '+peaceOutputFile+'.out '+dirName+'/')
			os.system('mv mstFile.mst '+dirName+'/')
			self.getPeaceResults(dirName, peaceOutputFile, peaceJobID)
		
		if self.runningWcd:
			os.system('mv wcd_'+estOutputFile+'.txt '+dirName+'/')
			if not self.wcdMatrixOut:
				# can't analyze results of matrix, only clusters
				self.getWcdResults(dirName, estOutputFile, wcdJobID)
		
		os.system('rm '+dirName+'/'+estOutputFile+'_dict.txt')
		
		print 'Script finished'

	def getPeaceResults(self, dirName, peaceOutputFile, jobID):
		oFile = open(dirName+'/'+peaceOutputFile+'.out', 'r')
		# Analyze and print results
		truePos = 0
		trueNeg = 0
		falsePos = 0
		falseNeg = 0
		currentCluster = {}
		geneCounts = {}
		geneTruePositives = {}
		estCount = 0
		
		for line in oFile:
			if line[0] == 'g':
				# an EST
				gTag = line.split('_')[0]
				estCount += 1
				val = currentCluster.get(gTag, 0)
				currentCluster[gTag] = val+1
				count = geneCounts.get(gTag, 0)
				geneCounts[gTag] = count+1
			elif line[0] == 'C':
				# new cluster
				# calc currentCluster, then reset
				iter = currentCluster.iteritems()
				thisCluTruePos = 0
				total = 0
				while True:
					try:
						cur = iter.next()
						oldVal = geneTruePositives.get(cur[0], 0)
						newVal = ((cur[1])*(cur[1]-1))/2
						geneTruePositives[cur[0]] = oldVal+newVal
						thisCluTruePos+=newVal
						total+=cur[1]
					except StopIteration:
						break
				truePos+=thisCluTruePos
				falsePos+=((total*(total-1))/2)-thisCluTruePos
				currentCluster = {}
		oFile.close()
		
		# calc currentCluster for the final time
		iter = currentCluster.iteritems()
		thisCluTruePos = 0
		total = 0
		while True:
			try:
				cur = iter.next()
				oldVal = geneTruePositives.get(cur[0], 0)
				newVal = ((cur[1])*(cur[1]-1))/2
				geneTruePositives[cur[0]] = oldVal+newVal
				thisCluTruePos+=newVal
				total+=cur[1]
			except StopIteration:
				break
			truePos+=thisCluTruePos
			falsePos+=((total*(total-1))/2)-thisCluTruePos
		# calc false negatives
		iter = geneCounts.iteritems()
		while True:
			try:
				cur = iter.next()
				falseNeg+=((cur[1]*(cur[1]-1))/2)-geneTruePositives.get(cur[0], 0)
			except StopIteration:
				break
				
		# calc true negatives (using n choose 2 minus other counts)
		trueNeg = ((estCount * (estCount-1))/2) - truePos - falsePos - falseNeg
		
		# calculate the time it took to run peace - note, take this out and replace with "time" call
		#runTime = os.stat(dirName+'/peace_script.o'+jobID)[8] - startTime

		# calculate Younden index (form: (a/(a+b) + d/(c+d) - 1))
		younden = -1.0
		if (truePos+falseNeg != 0):
			younden+=(truePos/(1.0*(truePos+falseNeg)))
		if (falsePos+trueNeg != 0):
			younden+=(trueNeg/(1.0*(falsePos+trueNeg)))

		rand = (1.0*(truePos+trueNeg))/(1.0*(falsePos+falseNeg+truePos+trueNeg))
	
		# Write analysis file to new directory
		aFile = open(dirName+'/analysis.txt', 'w')
		aFile.write('Peace Results\n')		
		aFile.write('True Positives: %d\n' %(truePos))
		aFile.write('True Negatives: %d\n' %(trueNeg))
		aFile.write('False Positives: %d\n' %(falsePos))
		aFile.write('False Negatives: %d\n' %(falseNeg))
		aFile.write('Younden Index: %f\n' %(younden))
		aFile.write('Uncorrected Rand Index: %f\n' %(rand))
		#aFile.write('Runtime: %f\n' %(runTime))
		aFile.write('\n')
		aFile.close()

		print 'Peace Results'
		print 'True Positives: %d' %(truePos)
		print 'True Negatives: %d' %(trueNeg)
		print 'False Positives: %d' %(falsePos)
		print 'False Negatives: %d' %(falseNeg)	
		print 'Younden Index: %f' %(younden)
		print 'Uncorrected Rand Index: %f' %(rand)
		#print 'Runtime: %f' %(runTime)
	# end getPeaceResults
	
	
	
	def getWcdResults(self, dirName, estOutputFile, jobID):		
		oFile = open(dirName+'/wcd_'+estOutputFile+'.txt', 'r')
		# Analyze and print results
		truePos = 0
		trueNeg = 0
		falsePos = 0
		falseNeg = 0
		currentCluster = {}
		geneCounts = {}
		geneTruePositives = {}
		estCount = 0
		
		for line in oFile:
			if line[0] == 'g':
				# an EST
				gTag = line.split('_')[0]
				estCount += 1
				val = currentCluster.get(gTag, 0)
				currentCluster[gTag] = val+1
				count = geneCounts.get(gTag, 0)
				geneCounts[gTag] = count+1
			elif line[0] == 'C':
				# new cluster
				# calc currentCluster, then reset
				iter = currentCluster.iteritems()
				thisCluTruePos = 0
				total = 0
				while True:
					try:
						cur = iter.next()
						oldVal = geneTruePositives.get(cur[0], 0)
						newVal = ((cur[1])*(cur[1]-1))/2
						geneTruePositives[cur[0]] = oldVal+newVal
						thisCluTruePos+=newVal
						total+=cur[1]
					except StopIteration:
						break
				truePos+=thisCluTruePos
				falsePos+=((total*(total-1))/2)-thisCluTruePos
				currentCluster = {}
		oFile.close()
		
		# calc currentCluster for the final time
		iter = currentCluster.iteritems()
		thisCluTruePos = 0
		total = 0
		while True:
			try:
				cur = iter.next()
				oldVal = geneTruePositives.get(cur[0], 0)
				newVal = ((cur[1])*(cur[1]-1))/2
				geneTruePositives[cur[0]] = oldVal+newVal
				thisCluTruePos+=newVal
				total+=cur[1]
			except StopIteration:
				break
			truePos+=thisCluTruePos
			falsePos+=((total*(total-1))/2)-thisCluTruePos
		# calc false negatives
		iter = geneCounts.iteritems()
		while True:
			try:
				cur = iter.next()
				falseNeg+=((cur[1]*(cur[1]-1))/2)-geneTruePositives.get(cur[0], 0)
			except StopIteration:
				break
				
		# calc true negatives (using n choose 2 minus other counts)
		trueNeg = ((estCount * (estCount-1))/2) - truePos - falsePos - falseNeg
		
		# calculate the time it took to run peace - note, take this out and replace with "time" call
		#runTime = os.stat(dirName+'/peace_script.o'+jobID)[8] - startTime

		# calculate Younden index (form: (a/(a+b) + d/(c+d) - 1))
		younden = -1.0
		if (truePos+falseNeg != 0):
			younden+=(truePos/(1.0*(truePos+falseNeg)))
		if (falsePos+trueNeg != 0):
			younden+=(trueNeg/(1.0*(falsePos+trueNeg)))

		rand = (1.0*(truePos+trueNeg))/(1.0*(falsePos+falseNeg+truePos+trueNeg))
	
		
		# Write analysis file to new directory
		aFile = open(dirName+'/analysis.txt', 'a')
		aFile.write('WCD Results\n')
		aFile.write('True Positives: %d\n' %(truePos))
		aFile.write('True Negatives: %d\n' %(trueNeg))
		aFile.write('False Positives: %d\n' %(falsePos))
		aFile.write('False Negatives: %d\n' %(falseNeg))
		aFile.write('Younden Index: %f\n' %(younden))
		aFile.write('Uncorrected Rand Index: %f\n' %(rand))
		#aFile.write('Runtime: %f\n' %(runTime))
		aFile.write('\n')
		aFile.close()
	
		print 'WCD Results'
		print 'True Positives: %d' %(truePos)
		print 'True Negatives: %d' %(trueNeg)
		print 'False Positives: %d' %(falsePos)
		print 'False Negatives: %d' %(falseNeg)
		print 'Younden Index: %f\n' %(younden)
		print 'Uncorrected Rand Index: %f\n' %(rand)
		#print 'Runtime: %f' %(runTime)
	# end getWCDResults

	def createJobFile(self, scriptName, proc, invocs, loadModules=False):
		if proc > 2:
			ppn = 2
			nodes = proc / 2
		else:
			ppn = proc
			nodes = 1

		# create the job file
		sFile = open(scriptName+'.job', 'w')
	
		sFile.write('#!/bin/bash -l\n')
		sFile.write('#PBS -N '+scriptName+'\n')
	
		sFile.write('#PBS -l nodes=%d:ppn=%d\n' %(nodes, ppn))
		sFile.write('#PBS -l walltime=100:00:00\n')
		sFile.write('#PBS -j oe\n')
	
		if loadModules:
			sFile.write('module load java\n')
			sFile.write('module load python\n')
		if self.runningPeace:
			sFile.write('export VIADEV_DEFAULT_RETRY_COUNT=7\n')
			sFile.write('export VIADEV_DEFAULT_TIME_OUT=21\n')
			#sFile.write('module load valgrind\n')

		# change to the correct directory
		sFile.write('cd '+os.getcwd()+'\n')

		for invoc in invocs:
			sFile.write(invoc+'\n')
			
		sFile.close()
	
		return scriptName+'.job'
	# end createJobFile

# start of main script
def main():
	# parse command line arguments

	usage= "usage: %prog [options] dataf cutf estsimparams"
	parser = OptionParser(usage=usage)
	parser.add_option("-p", "--peace", action="store_true", dest="runPeace", 
 			help="run peace with the default parameters.  Alternate parameters currently not supported.")
	parser.add_option("-w", "--wcd", 
			action="store_true", dest="runWcd",
			help="run wcd with the default parameters.  Alternate parameters currently not supported.")
	# TODO enable wcd parallel, and handle case where proc*2 > max 
	# procs allowed (when doing peace vs. wcd)
	parser.add_option("-c", "--processors", type="int", dest="proc",
	                  metavar=" PROC", help="request PROC processors")
	parser.add_option("-e", "--estfile", action="store_true", dest="specifyESTs",
			  help="specify an input file containing ESTs (don't create new simulated ESTs)")
	parser.set_defaults(runWcd=False, specifyESTs=False, runPeace=False, proc=1)
	(options, args) = parser.parse_args()
	if (len(args) < 2 and not options.specifyESTs) or (len(args) != 1 and options.specifyESTs):
	        parser.error("Incorrect number of arguments")

	# TODO there's no reason to do this.  request 1 extra and then 
	# invoke mpiexec -n proc
	proc = options.proc
	if proc <= 0:
		parser.error("Invalid argument, must have 1 or more processors")
	elif proc % 2 == 1 and proc != 1:
        	print "Odd number of processors specified, will be incremented by 1"
        	proc = proc+1

	if options.specifyESTs:
		pl = Pipeline(options.specifyESTs, args[0], None, 
			options.runPeace, options.runWcd)
	else:
		pl = Pipeline(options.specifyESTs, args[0], args[1], 
			options.runPeace, options.runWcd)
		pl.setESTParams(args[2:])
	pl.setProcessors(proc)
	pl.run()

if __name__ == "__main__":
	main()
