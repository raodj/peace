import os, sys, random

# Formatting script: takes any number of .fa files as arguments.
# Should work whether they are in ESTSim output format (i.e. full sequence on one line) or split up into multiple lines.
# Splits the sequences up, making them 80 characters to a line, and outputs them to new files.
# Also now randomly shuffles the ESTs and keeps track of their source genes (for the analysis routine later on).

# If peace format is used for wcd's output, the dictionary here is not needed.

for arg in sys.argv:
	newfn = arg[0:-3]
	ext = arg[-3:len(arg)]
	if ext == '.fa': # only run this on .fa files
		fa = open(arg, 'r')
		sqNum = -1
		sqs = dict()
		sqNums = []
		for line in fa:
			if line[0] == '>':
				sqNum+=1
				sqNums.append(sqNum)
				sqs[sqNum] = line+'\n'
			else:
				sqs[sqNum] = sqs[sqNum]+line.strip()
	
		random.shuffle(sqNums)
	
		output = open(newfn+'_fmt.fa', 'w')
		outputSqNum = 0
		outputDict = dict()

		for num in sqNums:
			line = sqs[num].split('\n')
			identifierLine = line[0]
			sequenceLine = line[2]

			# prefix it with a g (turns out the index by itself caused problems when using regex patterns later on)
			output.write('>g'+identifierLine[1:]+'\n')
			outputSqNumStr = '%d' %(outputSqNum)
			outputDict[outputSqNumStr] = 'g'+identifierLine[1:].split("_")[0]
			
			# write the sequence line
			i = 0
			while i+79 < len(sequenceLine):
				output.write(sequenceLine[i:i+79]+'\n')
				i = i + 79
			if i+1 < len(sequenceLine):
				output.write(sequenceLine[i:]+'\n')

			outputSqNum+=1

		output.close()
		fa.close()

		# write out the dictionary
		dictFile = open(newfn+'_dict.txt', 'w')
		dictFile.write(', '.join(outputDict.keys())+'\n')
		dictFile.write(', '.join(outputDict.values())+'\n')
		dictFile.close()
        	os.remove(arg) # clean up by removing the old .fa file
