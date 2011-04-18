import os, sys, time, popen2, random
from optparse import OptionParser



class Pipeline:
	def __init__(self, specifyESTs, faFile, digFile, peace, wcd, cap3):
		self.specifyESTs = specifyESTs
		self.runningPeace = peace
		self.runningWcd = wcd
		self.runningCap3 = cap3
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
			#peaceInvoc+='--estFile '+estOutputFile+'_fmt.fa --output '+peaceOutputFile+'.out --output-mst-file '+peaceOutputFile+'.mst --filters null'
			peaceInvoc+='--estFile '+estOutputFile+'_fmt.fa --output '+peaceOutputFile+'.out --output-mst-file '+peaceOutputFile+'.mst'

	
		if self.runningWcd:
			dirName = 'wcd_'+estOutputFile
			wcdInvoc = 'time mpiexec ./wcd -b -s'
			wcdInvoc += ' -o wcd_'+estOutputFile+'.txt '+estOutputFile+'_fmt.fa'

		if self.runningCap3:
			cap3OutputFile = estOutputFile+'.cap3'
			cap3Invoc = 'time cap3 '+estOutputFile+'_fmt.fa > '+cap3OutputFile
			# Assume we already have a directory for simplicity
	
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

		# first randomize script name
		# random auto-inits based on current system time
		scriptNum = random.randrange(0, 999999)


		if self.runningPeace and self.runningWcd:
		        # need to create and run 3 jobs
		        # estsim will run first
		        self.createJobFile('estsim_scr%d' %(scriptNum), 1, estSimInvocs, True)
		        p = popen2.popen2('qsub estsim_scr%d.job' %(scriptNum))
			jobID = p[0].read().split('.')[0]
			# Wait for script jobs to finish by checking for the output file
			while not os.path.exists('estsim_scr%d.o' %(scriptNum) +jobID):
		               	time.sleep(10)
		        # now create both peace and wcd job files, and run them in parallel
		        self.createJobFile('wcd_script%d' %(scriptNum), self.proc, [wcdInvoc])
		        self.createJobFile('peace_scr%d' %(scriptNum), self.proc, [peaceInvoc])
		        startTime = time.time()
			p = popen2.popen2('qsub peace_scr%d.job' %(scriptNum))
			peaceJobID = p[0].read().split('.')[0]
			p = popen2.popen2('qsub wcd_script%d.job' %(scriptNum))
			wcdJobID = p[0].read().split('.')[0]
			print "JobID for peace_script: "+peaceJobID
			print "JobID for wcd_script: "+wcdJobID
			
			# hacking this in
			if self.runningCap3:
				self.createJobFile('cap3_scr%d' %(scriptNum), 1, [cap3Invoc])
				p = popen2.popen2('qsub cap3_scr%d.job' %(scriptNum))
				cap3JobID = p[0].read().split('.')[0]
				print "JobID for cap3_script: "+cap3JobID

			while not os.path.exists('peace_scr%d.o' %(scriptNum) +peaceJobID) or not os.path.exists('wcd_script%d.o' %(scriptNum)+wcdJobID) or (self.runningCap3 and not os.path.exists('cap3_scr%d.o' %(scriptNum)+cap3JobID)):
				#if (time.time() - startTime) > 120:
				#	if not os.path.exists('peace_script.o'+peaceJobID):
				#		os.system('qdel '+peaceJobID)
				#		os.system('qsig -s 0 '+peaceJobID)
				#	if not os.path.exists('wcd_script.o'+wcdJobID):
				#		os.system('qdel '+wcdJobID)
				#		os.system('qsig -s 0 '+wcdJobID)
		                time.sleep(10)
		
		else:
			scriptName = 'estsim_scr%d' %(scriptNum)
		        if self.runningPeace:
		                estSimInvocs.append(peaceInvoc)
		                scriptName = 'peace_scr%d' %(scriptNum)
		        elif self.runningWcd:
		                estSimInvocs.append(wcdInvoc)
		                scriptName = 'wcd_script%d' %(scriptNum)
		        self.createJobFile(scriptName, self.proc, estSimInvocs, True)
		        startTime = time.time()
			if self.proc == 1:
			        p = popen2.popen2('qsub -q serial_restr '+scriptName+'.job')
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
			#os.system('mv '+estOutputFile+'_dict.txt '+dirName+'/')
		
		if self.runningPeace and self.runningWcd:
			os.system('mv estsim_scr%d.o' %(scriptNum)+jobID+' '+dirName+'/')
			os.system('mv peace_scr%d.o' %(scriptNum) +peaceJobID+' '+dirName+'/')
			os.system('mv wcd_script%d.o' %(scriptNum) +wcdJobID+' '+dirName+'/')
			os.system('rm peace_scr%d.job' %(scriptNum))
			os.system('rm estsim_scr%d.job' %(scriptNum))
			os.system('rm wcd_script%d.job' %(scriptNum))
		else:
			#os.system('rm '+scriptName+'.job')
			os.system('mv '+scriptName+'.o'+jobID+' '+dirName+'/')
			peaceJobID = jobID
			wcdJobID = jobID

		if self.runningPeace:
			os.system('mv '+peaceOutputFile+'.out '+dirName+'/')
			os.system('mv '+peaceOutputFile+'.mst '+dirName+'/')
			self.getPeaceResults(dirName, peaceOutputFile, peaceJobID, scriptNum)
		
		if self.runningWcd:
			os.system('mv wcd_'+estOutputFile+'.txt '+dirName+'/')
			if not self.wcdMatrixOut:
				# can't analyze results of matrix, only clusters
				self.getWcdResults(dirName, estOutputFile, wcdJobID, scriptNum)

		if self.runningCap3:
			os.system('mv cap3_scr%d.o' %(scriptNum) + cap3JobID+' '+dirName+'/')
			os.system('rm cap3_scr%d.job' %(scriptNum))
			os.system('mv '+cap3OutputFile+' '+dirName+'/')
			cap3SingletFile = estOutputFile+'_fmt.fa.cap.singlets'
			os.system('mv '+cap3SingletFile+' '+dirName+'/')
			os.system('rm '+estOutputFile+'_fmt.fa.cap.ace')
			os.system('rm '+estOutputFile+'_fmt.fa.cap.contigs')
			os.system('rm '+estOutputFile+'_fmt.fa.cap.contigs.links')
			os.system('rm '+estOutputFile+'_fmt.fa.cap.contigs.qual')
			os.system('rm '+estOutputFile+'_fmt.fa.cap.info')

			self.getCap3Results(dirName, cap3OutputFile, cap3SingletFile, cap3JobID, scriptNum)		
		
		#os.system('rm '+dirName+'/'+estOutputFile+'_dict.txt')
		
		print 'Script finished'

	def getPeaceResults(self, dirName, peaceOutputFile, jobID, num):
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
		
		# get runtime
		runTime = 'n/a'
		rtFile = open(dirName+'/peace_scr%d.o' %(num)+jobID, 'r')
		for line in rtFile:
			rtList = line.split()
			if len(rtList) >= 2 and rtList[0] == 'real':
				runTime = rtList[1]
		rtFile.close()

		# calculate Younden index (form: (a/(a+b) + d/(c+d) - 1))
		younden = -1.0
		if (truePos+falseNeg != 0):
			younden+=(truePos/(1.0*(truePos+falseNeg)))
		if (falsePos+trueNeg != 0):
			younden+=(trueNeg/(1.0*(falsePos+trueNeg)))

		rand = (1.0*(truePos+trueNeg))/(1.0*(falsePos+falseNeg+truePos+trueNeg))

		sensitivity = (1.0*truePos)/(1.0*(truePos+falseNeg))

		jaccard = (1.0*truePos)/(1.0*(truePos+falseNeg+falsePos))
	
		# Write analysis file to new directory
		aFile = open(dirName+'/analysis.txt', 'w')
		aFile.write('Peace Results\n')		
		aFile.write('True Positives: %d\n' %(truePos))
		aFile.write('True Negatives: %d\n' %(trueNeg))
		aFile.write('False Positives: %d\n' %(falsePos))
		aFile.write('False Negatives: %d\n' %(falseNeg))
		aFile.write('Younden Index: %f\n' %(younden))
		aFile.write('Uncorrected Rand Index: %f\n' %(rand))
		aFile.write('Sensitivity: %f\n' %(sensitivity))
		aFile.write('Jaccard Index: %f\n' %(jaccard))
		aFile.write('Runtime: '+runTime+'\n')
		aFile.write('\n')
		aFile.close()

		print 'Peace Results'
		print 'True Positives: %d' %(truePos)
		print 'True Negatives: %d' %(trueNeg)
		print 'False Positives: %d' %(falsePos)
		print 'False Negatives: %d' %(falseNeg)	
		print 'Younden Index: %f' %(younden)
		print 'Uncorrected Rand Index: %f' %(rand)
		print 'Sensitivity: %f' %(sensitivity)
		print 'Jaccard Index: %f' %(jaccard)
		print 'Runtime: '+runTime
	# end getPeaceResults
	
	
	
	def getWcdResults(self, dirName, estOutputFile, jobID, num):
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
		
		# get runtime
		runTime = 'n/a'
		rtFile = open(dirName+'/wcd_script%d.o' %(num)+jobID, 'r')
                for line in rtFile:
                        rtList = line.split()
                        if len(rtList) >= 2 and rtList[0] == 'real':
                                runTime = rtList[1]
                rtFile.close()

		# calculate Younden index (form: (a/(a+b) + d/(c+d) - 1))
		younden = -1.0
		if (truePos+falseNeg != 0):
			younden+=(truePos/(1.0*(truePos+falseNeg)))
		if (falsePos+trueNeg != 0):
			younden+=(trueNeg/(1.0*(falsePos+trueNeg)))

		rand = (1.0*(truePos+trueNeg))/(1.0*(falsePos+falseNeg+truePos+trueNeg))
	
		sensitivity = (1.0*truePos)/(1.0*(truePos+falseNeg))

		jaccard = (1.0*truePos)/(1.0*(truePos+falseNeg+falsePos))
		
		# Write analysis file to new directory
		aFile = open(dirName+'/analysis.txt', 'a')
		aFile.write('WCD Results\n')
		aFile.write('True Positives: %d\n' %(truePos))
		aFile.write('True Negatives: %d\n' %(trueNeg))
		aFile.write('False Positives: %d\n' %(falsePos))
		aFile.write('False Negatives: %d\n' %(falseNeg))
		aFile.write('Younden Index: %f\n' %(younden))
		aFile.write('Uncorrected Rand Index: %f\n' %(rand))
		aFile.write('Sensitivity: %f\n' %(sensitivity))
		aFile.write('Jaccard Index: %f\n' %(jaccard))
		aFile.write('Runtime: '+runTime+'\n')
		aFile.write('\n')
		aFile.close()
	
		print 'WCD Results'
		print 'True Positives: %d' %(truePos)
		print 'True Negatives: %d' %(trueNeg)
		print 'False Positives: %d' %(falsePos)
		print 'False Negatives: %d' %(falseNeg)
		print 'Younden Index: %f\n' %(younden)
		print 'Uncorrected Rand Index: %f\n' %(rand)
		print 'Sensitivity: %f' %(sensitivity)
		print 'Jaccard Index: %f' %(jaccard)
		print 'Runtime: '+runTime
	# end getWCDResults

	def getCap3Results(self, dirName, cap3OutputFile, cap3SingletFile, jobID, num):
		oFile = open(dirName+'/'+cap3OutputFile, 'r')
		singFile = open(dirName+'/'+cap3SingletFile, 'r')
		# Analyze and print results
		truePos = 0
		trueNeg = 0
		falsePos = 0
		falseNeg = 0
		currentCluster = {}
		geneCounts = {}
		geneTruePositives = {}
		estCount = 0

		# Do singletons first
		for line in singFile:
			if line[0] == '>':
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
				# an EST
				gTag = line.split('_')[0]
				estCount += 1
				val = currentCluster.get(gTag, 0)
				currentCluster[gTag] = val+1
				count = geneCounts.get(gTag, 0)
				geneCounts[gTag] = count+1
			# else, do nothing

		# final currentCluster will be calculated for first cluster in contig file
		singFile.close()

		# do contig file
		for line in oFile:
			if line.strip() == "DETAILED DISPLAY OF CONTIGS":
				# Got all the clusters
				break
			elif len(line.strip()) > 0 and line.strip()[0] == 'g':
				# an EST
				gTag = line.strip().split('_')[0]
				estCount += 1
				val = currentCluster.get(gTag, 0)
				currentCluster[gTag] = val+1
				count = geneCounts.get(gTag, 0)
				geneCounts[gTag] = count+1
			elif len(line) > 0 and line[0] == '*':
				# new contig
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

		# get runtime
		runTime = 'n/a'
		rtFile = open(dirName+'/cap3_scr%d.o' %(num)+jobID, 'r')
		for line in rtFile:
			rtList = line.split()
			if len(rtList) >= 2 and rtList[0] == 'real':
				runTime = rtList[1]
		rtFile.close()

		# calculate Younden index (form: (a/(a+b) + d/(c+d) - 1))
		younden = -1.0
		if (truePos+falseNeg != 0):
			younden+=(truePos/(1.0*(truePos+falseNeg)))
		if (falsePos+trueNeg != 0):
			younden+=(trueNeg/(1.0*(falsePos+trueNeg)))

		rand = (1.0*(truePos+trueNeg))/(1.0*(falsePos+falseNeg+truePos+trueNeg))

		sensitivity = (1.0*truePos)/(1.0*(truePos+falseNeg))

		jaccard = (1.0*truePos)/(1.0*(truePos+falseNeg+falsePos))
		
		# Write analysis file to new directory
		aFile = open(dirName+'/analysis.txt', 'a')
		aFile.write('CAP3 Results\n')		
		aFile.write('True Positives: %d\n' %(truePos))
		aFile.write('True Negatives: %d\n' %(trueNeg))
		aFile.write('False Positives: %d\n' %(falsePos))
		aFile.write('False Negatives: %d\n' %(falseNeg))
		aFile.write('Younden Index: %f\n' %(younden))
		aFile.write('Uncorrected Rand Index: %f\n' %(rand))
		aFile.write('Sensitivity: %f\n' %(sensitivity))
		aFile.write('Jaccard Index: %f\n' %(jaccard))
		aFile.write('Runtime: '+runTime+'\n')
		aFile.write('\n')
		aFile.close()

		print 'CAP3 Results'
		print 'True Positives: %d' %(truePos)
		print 'True Negatives: %d' %(trueNeg)
		print 'False Positives: %d' %(falsePos)
		print 'False Negatives: %d' %(falseNeg)	
		print 'Younden Index: %f' %(younden)
		print 'Uncorrected Rand Index: %f' %(rand)
		print 'Sensitivity: %f' %(sensitivity)
		print 'Jaccard Index: %f' %(jaccard)
		print 'Runtime: '+runTime
	# end getCap3Results

	

	def createJobFile(self, scriptName, proc, invocs, loadModules=False):
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
		sFile.write('#PBS -l walltime=1:00:00\n')
		sFile.write('#PBS -j oe\n')
		if nodes == 1 and ppn == 1:
			sFile.write('#PBS -q serial_restr\n')

	
		if loadModules:
			sFile.write('module load java\n')
		if self.runningPeace:
			sFile.write('export VIADEV_DEFAULT_RETRY_COUNT=7\n')
			sFile.write('export VIADEV_DEFAULT_TIME_OUT=21\n')
			#sFile.write('module load valgrind\n')
		if self.runningCap3:
			sFile.write('module load cap3\n')

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
	print "--------------------- pipeline.py running -------------------\n"
	usage= "usage: %prog [options] dataf cutf estsimparams"
	parser = OptionParser(usage=usage)
	parser.add_option("-p", "--peace", action="store_true", dest="runPeace", 
 			help="run peace with the default parameters.  Alternate parameters currently not supported.")
	parser.add_option("-w", "--wcd", 
			action="store_true", dest="runWcd",
			help="run wcd with the default parameters.  Alternate parameters currently not supported.")
	parser.add_option("-d", "--cap3", action="store_true", dest="runCap3")
	parser.add_option("-c", "--processors", type="int", dest="proc",
	                  metavar=" PROC", help="request PROC processors")
	parser.add_option("-e", "--estfile", action="store_true", dest="specifyESTs",
			  help="specify an input file containing ESTs (don't create new simulated ESTs)")
	parser.set_defaults(runWcd=False, specifyESTs=False, runPeace=False, proc=1, runCap3=False)
	(options, args) = parser.parse_args()
	if (len(args) < 2 and not options.specifyESTs) or (len(args) != 1 and options.specifyESTs):
	        parser.error("Incorrect number of arguments")

	proc = options.proc
	if proc <= 0:
		parser.error("Invalid argument, must have 1 or more processors")
	#elif proc % 2 == 1 and proc != 1:
        #	print "Odd number of processors specified, will be 
	#incremented by 1"
        #	proc = proc+1

	if options.specifyESTs:
		pl = Pipeline(options.specifyESTs, args[0], None, 
			options.runPeace, options.runWcd, options.runCap3)
	else:
		pl = Pipeline(options.specifyESTs, args[0], args[1], 
			options.runPeace, options.runWcd, options.runCap3)
		pl.setESTParams(args[2:])
	pl.setProcessors(proc)
	pl.run()

if __name__ == "__main__":
	main()
