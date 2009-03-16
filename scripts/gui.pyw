# GUI for handling pipeline script input/output

from Tkinter import *
from Pipeline import Pipeline
import tkFileDialog
import tkMessageBox

class MainWindow(Frame):
	def __init__(self, parent):
		frame = Frame(parent)
		#super(MainWindow, self).__init__(parent)
		self.parent = parent
		frame.grid(row=0, column=0)
		
		# variables
		
		self.faFN=StringVar()
		self.faFN.set("")
		self.digFN=StringVar()
		self.digFN.set("")
		self.runPeace=IntVar()
		self.runWcd=IntVar()
		self.runESTSim=IntVar()
		self.runESTSim.set(1)
		self.wcdMatrixOutput = IntVar()
		self.analyzer=StringVar()
		self.analyzer.set("clu")
		self.clusterMaker=StringVar()
		self.clusterMaker.set("mst")
		self.estIdx=IntVar()
		self.estIdx.set(0)
		self.procs=IntVar()
		self.procs.set(1)
		self.analyzerB = []
		self.clusterMakerB = []
		self.numClones=IntVar()
		self.numClones.set(1)

		self.sbe=IntVar()
		self.sbe.set(0)
		self.stutter=IntVar()
		self.stutter.set(0)
		self.ligate=IntVar()
		self.ligate.set(0)
		self.alpha=DoubleVar()
		self.alpha.set(0)
		self.beta=IntVar()
		self.beta.set(0)
		self.gamma=DoubleVar()
		self.gamma.set(0)
		self.zeta=DoubleVar()
		self.zeta.set(0)
		self.xi=DoubleVar()
		self.xi.set(0)
		self.kappa=IntVar()
		self.kappa.set(0)
		self.lambda_=IntVar()
		self.lambda_.set(0)
		self.mu=IntVar()
		self.mu.set(0)
		self.nu=IntVar()
		self.nu.set(0)
		self.stuttering=DoubleVar()
		self.stuttering.set(0)
		self.ligations=DoubleVar()
		self.ligations.set(0)
		
		# objects in form
		
		# beginning buttons/labels
		
		runEstSimButton = Checkbutton(frame, text="Generate Simulated ESTs?", variable=self.runESTSim, command=self.toggleESTSim)
		self.faFileButton = Button(frame,text="Sequence Filename", anchor=W, command=self.getfaFN)
		faFileLabel = Label(frame, textvariable=self.faFN, relief=SUNKEN, anchor=E)
		self.digFileButton = Button(frame,text="Cut Filename", anchor=W, command=self.getdigFN)
		digFileLabel = Label(frame, textvariable=self.digFN, relief=SUNKEN, anchor=E)

		# peace and peace parameters
		
		runPeaceButton = Checkbutton(frame, text="Run Peace?", variable=self.runPeace, command=self.togglePeace)
		analyzerLabel = Label(frame, text="Analyzer:", anchor=W)
		clusterLabel = Label(frame, text="Cluster Maker:", anchor=W)
		estIdxLabel = Label(frame, text="estIdx:", anchor=W)
		
		self.analyzerB.append(Radiobutton(frame, text="clu", variable=self.analyzer, value="clu", takefocus=FALSE, state=DISABLED))
		self.clusterMakerB.append(Radiobutton(frame, text="mst", variable=self.clusterMaker, value="mst", takefocus=FALSE, state=DISABLED))
	
		self.estIdxEntry = Entry(frame, textvariable=self.estIdx, takefocus=FALSE, state=DISABLED)

		procLabel = Label(frame, text="Processors:", anchor=W)
		self.procEntry = Entry(frame, textvariable=self.procs, takefocus=FALSE, state=DISABLED)
		
		estSimParamsLabel = Label(frame, text="ESTSim Parameters", relief=RIDGE, anchor=CENTER)

		otherParamsLabel = Label(frame, text="Peace and WCD Parameters", relief=RIDGE, anchor=CENTER)
		
		# wcd and wcd parameters
		
		runWcdButton = Checkbutton(frame, text="Run WCD?", variable=self.runWcd, command=self.toggleWcd)
		self.wcdMatrixButton = Checkbutton(frame, text="WCD Matrix Output?", variable=self.wcdMatrixOutput, takefocus=FALSE, state=DISABLED)
		
		# estsim parameters

		numClonesLabel = Label(frame, text="Number of EST Clones:", anchor=W)
		self.numClonesEntry = Entry(frame, textvariable=self.numClones)

		self.sbeButton = Checkbutton(frame, text="Use sbe fault model?", variable=self.sbe, command=self.toggleSbe)
		self.alphaEntry = Entry(frame, textvariable=self.alpha, takefocus=FALSE, state=DISABLED)
		self.betaEntry = Entry(frame, textvariable=self.beta, takefocus=FALSE, state=DISABLED)
		self.gammaEntry = Entry(frame, textvariable=self.gamma, takefocus=FALSE, state=DISABLED)
		self.zetaEntry = Entry(frame, textvariable=self.zeta, takefocus=FALSE, state=DISABLED)
		self.xiEntry = Entry(frame, textvariable=self.xi, takefocus=FALSE, state=DISABLED)
		self.kappaEntry = Entry(frame, textvariable=self.kappa, takefocus=FALSE, state=DISABLED)
		self.lambdaEntry = Entry(frame, textvariable=self.lambda_, takefocus=FALSE, state=DISABLED)
		self.muEntry = Entry(frame, textvariable=self.mu, takefocus=FALSE, state=DISABLED)
		self.nuEntry = Entry(frame, textvariable=self.nu, takefocus=FALSE, state=DISABLED)
		alphaLabel = Label(frame, text="Single-Base Error Rate:", anchor=W)
		betaLabel = Label(frame, text="Primer Interference Margin:", anchor=W)
		gammaLabel = Label(frame, text="Primer Interference Rate:", anchor=W)
		zetaLabel = Label(frame, text="Rapid Polymerase Decay Rate:", anchor=W)
		xiLabel = Label(frame, text="Gentle Polymerase Decay Rate:", anchor=W)
		kappaLabel = Label(frame, text="Proportion of Substitutions:", anchor=W)
		lambdaLabel = Label(frame, text="Proportion of Deletions:", anchor=W)
		muLabel = Label(frame, text="Proportion of Insertions:", anchor=W)
		nuLabel = Label(frame, text="Proportion of Ns:", anchor=W)
		self.stutterButton = Checkbutton(frame, text="Use stutter fault model?", variable=self.stutter, command=self.toggleStutter)
		stutterLabel = Label(frame, text="Stuttering:", anchor=W)
		self.stutterEntry = Entry(frame, textvariable=self.stuttering, takefocus=FALSE, state=DISABLED)
		self.ligateButton = Checkbutton(frame, text="Use ligate fault model?", variable=self.ligate, command=self.toggleLigate)
		ligateLabel = Label(frame, text="Ligation:", anchor=W)
		self.ligateEntry = Entry(frame, textvariable=self.ligations, takefocus=FALSE, state=DISABLED)
		
		test = Button(frame, text="?", anchor=E, takefocus=FALSE, command=self.sbeHelp)

		self.startB = Button(frame, text="Run Pipeline Script", fg="red", anchor=CENTER, command=self.start, takefocus=FALSE, state=DISABLED)
		
		# grid setup
		
		runEstSimButton.grid(row=0, column=0, columnspan=2, padx=2, pady=2, sticky=EW)
		
		self.faFileButton.grid(row=1, column=0, padx=2, pady=2, sticky=W)
		faFileLabel.grid(row=1, column=1, padx=2, pady=2, sticky=EW)
		
		self.digFileButton.grid(row=2, column=0, padx=2, pady=2, sticky=W)
		digFileLabel.grid(row=2, column=1, padx=2, pady=2, sticky=EW)
		
		otherParamsLabel.grid(row=3, column=0, columnspan=2, padx=2, pady=2, sticky=EW)

		runPeaceButton.grid(row=4, column=0, columnspan=2, padx=2, pady=2, sticky=EW)
		
		analyzerLabel.grid(row=5, column=0, padx=2, pady=2, sticky=W)
		i = 1
		for button in self.analyzerB:
			button.grid(row=5, column=i, padx=2, pady=2, sticky=EW)
			i+=1
		
		i = 1
		clusterLabel.grid(row=6, column=0, padx=2, pady=2, sticky=W)
		for button in self.clusterMakerB:
			button.grid(row=6, column=i, padx=2, pady=2, sticky=EW)
			i+=1

		estIdxLabel.grid(row=7, column=0, padx=2, pady=2, sticky=W)
		self.estIdxEntry.grid(row=7, column=1, padx=2, pady=2, sticky=EW)
			
		procLabel.grid(row=8, column=0, padx=2, pady=2, sticky=W)
		self.procEntry.grid(row=8, column=1, padx=2, pady=2, sticky=EW)
		
		runWcdButton.grid(row=9, column=0, columnspan=2, padx=2, pady=2, sticky=EW)
		self.wcdMatrixButton.grid(row=10, column=0, padx=2, pady=2, sticky=EW)
		
		estSimParamsLabel.grid(row=11, column=0, columnspan=2, padx=2, pady=2, sticky=EW)

		numClonesLabel.grid(row=12, column=0, padx=2, pady=2, sticky=W)
		self.numClonesEntry.grid(row=12, column=1, padx=2, pady=2, sticky=EW)

		self.sbeButton.grid(row=13, column=0, columnspan=2, padx=2, pady=2, sticky=EW)
		self.alphaEntry.grid(row=14, column=1, padx=2, pady=2, sticky=EW)
		self.betaEntry.grid(row=15, column=1, padx=2, pady=2, sticky=EW)
		self.gammaEntry.grid(row=16, column=1, padx=2, pady=2, sticky=EW)
		self.zetaEntry.grid(row=17, column=1, padx=2, pady=2, sticky=EW)
		self.xiEntry.grid(row=18, column=1, padx=2, pady=2, sticky=EW)
		self.kappaEntry.grid(row=19, column=1, padx=2, pady=2, sticky=EW)
		self.lambdaEntry.grid(row=20, column=1, padx=2, pady=2, sticky=EW)
		self.muEntry.grid(row=21, column=1, padx=2, pady=2, sticky=EW)
		self.nuEntry.grid(row=22, column=1, padx=2, pady=2, sticky=EW)
		alphaLabel.grid(row=14, column=0, padx=2, pady=2, sticky=W)
		betaLabel.grid(row=15, column=0, padx=2, pady=2, sticky=W)
		gammaLabel.grid(row=16, column=0, padx=2, pady=2, sticky=W)
		zetaLabel.grid(row=17, column=0, padx=2, pady=2, sticky=W)
		xiLabel.grid(row=18, column=0, padx=2, pady=2, sticky=W)
		kappaLabel.grid(row=19, column=0, padx=2, pady=2, sticky=W)
		lambdaLabel.grid(row=20, column=0, padx=2, pady=2, sticky=W)
		muLabel.grid(row=21, column=0, padx=2, pady=2, sticky=W)
		nuLabel.grid(row=22, column=0, padx=2, pady=2, sticky=W)

		self.stutterButton.grid(row=23, column=0, columnspan=2, padx=2, pady=2, sticky=EW)
		self.stutterEntry.grid(row=24, column=1, padx=2, pady=2, sticky=EW)
		stutterLabel.grid(row=24, column=0, padx=2, pady=2, sticky=W)

		self.ligateButton.grid(row=25, column=0, columnspan=2, padx=2, pady=2, sticky=EW)
		self.ligateEntry.grid(row=26, column=1, padx=2, pady=2, sticky=EW)
		ligateLabel.grid(row=26, column=0, padx=2, pady=2, sticky=W)
		
		self.startB.grid(row=27, column=0, columnspan=2, padx=2, pady=2, sticky=EW)

		test.grid(row=13, column=2, padx=2, pady=2, sticky=E)

		parent.bind("<F1>", self.sbeHelp)
		#parent.bind("<Control-q>", self.quit)
		#parent.bind("<Escape>", self.quit)
		
	# functions
		
	def getfaFN(self):
		self.faFN.set(tkFileDialog.askopenfilename(filetypes=["{FastA Files} {.fasta}", "{FastA Files} {.fa}"], title="Select Sequences File"))
		self.enableStart()
		
	def getdigFN(self):
		self.digFN.set(tkFileDialog.askopenfilename(filetypes=["{Cut Files} {.dig}"], title="Select Cut File"))
		self.enableStart()
		
	def getestFN(self):
		self.faFN.set(tkFileDialog.askopenfilename(filetypes=["{FastA Files} {.fasta}", "{FastA Files} {.fa}"], title="Select EST File"))
		self.enableStart()
		
	def togglePeace(self):
		if self.runPeace.get() == 1:
			for button in self.analyzerB:
				button.config(state=NORMAL, takefocus=TRUE)
			for button in self.clusterMakerB:
				button.config(state=NORMAL, takefocus=TRUE)
			self.estIdxEntry.config(state=NORMAL, takefocus=TRUE)
			self.procEntry.config(state=NORMAL, takefocus=TRUE)
		else:
			for button in self.analyzerB:
				button.config(state=DISABLED, takefocus=FALSE)
			for button in self.clusterMakerB:
				button.config(state=DISABLED, takefocus=FALSE)
			self.estIdxEntry.config(state=DISABLED, takefocus=FALSE)
			self.procEntry.config(state=DISABLED, takefocus=FALSE)
			
	def toggleWcd(self):
		if self.runWcd.get() == 1:
			self.wcdMatrixButton.config(state=NORMAL, takefocus=TRUE)
		else:
			self.wcdMatrixButton.config(state=DISABLED, takefocus=FALSE)
				
	def toggleESTSim(self):
		# doesn't include estsim params
		if self.runESTSim.get() == 1:
			self.digFileButton.config(state=NORMAL, takefocus=TRUE)
			self.faFileButton.config(text="Sequence Filename", command=self.getfaFN)
			self.sbeButton.config(state=NORMAL, takefocus=TRUE)
			self.stutterButton.config(state=NORMAL, takefocus=TRUE)
			self.ligateButton.config(state=NORMAL, takefocus=TRUE)
			self.numClonesEntry.config(state=NORMAL, takefocus=TRUE)
		else:
			self.digFileButton.config(state=DISABLED, takefocus=FALSE)
			self.faFileButton.config(text="EST Filename", command=self.getestFN)
			self.sbeButton.config(state=DISABLED, takefocus=FALSE)
			self.stutterButton.config(state=DISABLED, takefocus=FALSE)
			self.ligateButton.config(state=DISABLED, takefocus=FALSE)
			self.numClonesEntry.config(state=DISABLED, takefocus=FALSE)
		# refresh these
		self.toggleSbe()
		self.toggleStutter()
		self.toggleLigate()
		self.enableStart()

	def toggleSbe(self):
		if self.sbe.get() == 1 and self.runESTSim.get() == 1:
			self.alphaEntry.config(state=NORMAL, takefocus=TRUE)
			self.betaEntry.config(state=NORMAL, takefocus=TRUE)
			self.gammaEntry.config(state=NORMAL, takefocus=TRUE)
			self.zetaEntry.config(state=NORMAL, takefocus=TRUE)
			self.xiEntry.config(state=NORMAL, takefocus=TRUE)
			self.kappaEntry.config(state=NORMAL, takefocus=TRUE)
			self.lambdaEntry.config(state=NORMAL, takefocus=TRUE)
			self.muEntry.config(state=NORMAL, takefocus=TRUE)
			self.nuEntry.config(state=NORMAL, takefocus=TRUE)
		else:
			self.alphaEntry.config(state=DISABLED, takefocus=FALSE)
			self.betaEntry.config(state=DISABLED, takefocus=FALSE)
			self.gammaEntry.config(state=DISABLED, takefocus=FALSE)
			self.zetaEntry.config(state=DISABLED, takefocus=FALSE)
			self.xiEntry.config(state=DISABLED, takefocus=FALSE)
			self.kappaEntry.config(state=DISABLED, takefocus=FALSE)
			self.lambdaEntry.config(state=DISABLED, takefocus=FALSE)
			self.muEntry.config(state=DISABLED, takefocus=FALSE)
			self.nuEntry.config(state=DISABLED, takefocus=FALSE)

	def toggleStutter(self):
		if self.stutter.get() == 1 and self.runESTSim.get() == 1:
			self.stutterEntry.config(state=NORMAL, takefocus=TRUE)
		else:
			self.stutterEntry.config(state=DISABLED, takefocus=FALSE)

	def toggleLigate(self):
		if self.ligate.get() == 1 and self.runESTSim.get() == 1:
			self.ligateEntry.config(state=NORMAL, takefocus=TRUE)
		else:
			self.ligateEntry.config(state=DISABLED, takefocus=FALSE)
		
	def enableStart(self):
		if self.runESTSim.get() == 1:
			if self.faFN.get() != "" and self.digFN.get() != "":
				self.startB.config(state=NORMAL, takefocus=TRUE)
			else:
				self.startB.config(state=DISABLED, takefocus=FALSE)
		else:
			if self.faFN.get() != "":
				self.startB.config(state=NORMAL, takefocus=TRUE)
			else:
				self.startB.config(state=DISABLED, takefocus=FALSE)
			
	# At present does not do error checking, just lets them put in whatever
	def start(self):
		# if error checking fails just display an error msg and return
		# otherwise this will tell the pipeline to run, and close itself
	
		if self.runESTSim.get() == 1:			
			pl = Pipeline(True, self.faFN, self.digFN)
			estParams = []
			estParams.append(self.numClones.get())
			# get all params from variables
			# error check all vars to be within bounds
			if self.sbe.get() == 1:
				estParams.append('sbe')
				estParams.append(self.alpha.get())
				estParams.append(self.beta.get())
				estParams.append(self.gamma.get())
				estParams.append(self.zeta.get())
				estParams.append(self.xi.get())
				estParams.append(self.kappa.get())
				estParams.append(self.lambda_.get())
				estParams.append(self.mu.get())
				estParams.append(self.nu.get())
			if self.stutter.get() == 1:
				estParams.append('stutter')
				estParams.append(self.stuttering.get())
			if self.ligate.get() == 1:
				estParams.append('ligate')
				estParams.append(self.ligations.get())
			# set params
			pl.setESTParams(estParams)
		else:
			pl = Pipeline(False, self.faFN, None)
			
		if self.runPeace.get() == 1:
			# must error check procs (see pipeline.py) and estIdx (can't be < 0)
			pl.setPeaceParams(self.analyzer.get(), self.clusterMaker.get(), self.estIdx.get(), self.proc.get())
		
		if self.runWcd.get() == 1:
			if self.wcdMatrixOutput.get() == 1:
				pl.setWcdParams(True)
			else:
				pl.setWcdParams(False)
		
		pl.run()
		
		# exit
		self.parent.destroy()
		
	def quit(self, event=None):
		self.parent.destroy()

	def sbeHelp(self, event=None):
		tkMessageBox.showinfo("sbe Help", "Models the single-base error rate.  Parameters:\n"+
		"Single-base error rate: overall rate of single-base errors (between 0 and 1).\n"+
		"Primer interference margin: margin at the beginning of the EST where errors should be more likely (set to 0 for none).\n"+
		"Primer interference rate: rate of error in the primer interference margin (set to 0 for no change).\n"+
		"Rapid polymerase decay rate: rate of error at the end of the EST (set to 0 for no change).\n"+
		"Gentle polymerase decay rate: rate of error increase as we approach the end of the EST (set to 0 for no change).\n"+
		"Proportions give the proportion of single-base errors for the given type of error.")
	
application = Tk()
application.title("Peace/WCD Pipeline Input")
window=MainWindow(application)
application.protocol("WM_DELETE_WINDOW",window.quit)
application.mainloop()