from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files

import os
from csv import reader

def SingleQubitRB(qubit, seqFile, showPlot=False):
	"""

	Single qubit randomized benchmarking using 90 and 180 generators. 

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel) 
	seqFile : file containing sequence strings
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""	

	#Setup a pulse library
	pulseLib = {'X90p':X90(qubit), 'X90m':X90m(qubit), 'Xp':X(qubit), 'Xm':Xm(qubit),
	            'Y90p':Y90(qubit), 'Y90m':Y90m(qubit), 'Yp':Y(qubit), 'Ym':Ym(qubit),
	            'QId':Id(qubit)}
	measBlock = MEAS(qubit)

	with open(seqFile,'r') as FID:
		fileReader = reader(FID, delimiter='\t')
		seqs = []
		for pulseSeqStr in fileReader:
			seq = []
			for pulseStr in pulseSeqStr:
				seq.append(pulseLib[pulseStr])
			seq.append(measBlock)
			seqs.append(seq)

	#Tack on the calibration scalings
	seqs += [[Id(qubit), measBlock], [Id(qubit), measBlock], [X(qubit), measBlock], [X(qubit), measBlock]]

	fileNames = compile_to_hardware(seqs, 'RB/RB')
	print(fileNames)

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin


def SingleQubitRB_AC(qubit, seqFile, showPlot=False):
	"""

	Single qubit randomized benchmarking using atomic Clifford pulses. 

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel) 
	seqFile : file containing sequence strings
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""	
	#Setup a pulse library
	pulseLib = [AC(qubit, cliffNum) for cliffNum in range(24)]
	pulseLib.append(pulseLib[0])
	measBlock = MEAS(qubit)

	with open(seqFile,'r') as FID:
		fileReader = reader(FID)
		seqs = []
		for pulseSeqStr in fileReader:
			seq = []
			for pulseStr in pulseSeqStr:
				seq.append(pulseLib[int(pulseStr)-1])
			seq.append(measBlock)
			seqs.append(seq)

	# #Tack on the calibration scalings
	# seqs += [[Id(qubit), measBlock], [Id(qubit), measBlock], [X(qubit), measBlock], [X(qubit), measBlock]]
	# print('Number of sequences: {0}'.format(len(seqs)))

	# fileNames = compile_to_hardware(seqs, 'RB/RB')
	# print fileNames
	#Hack for limited APS waveform memory and break it up into multiple files
	#We've shuffled the sequences so that we loop through each gate length on the inner loop
	numRandomizations = 36
	numGateLengths = 17
	for ct in range(numRandomizations):
		# chunk = seqs[2*numGateLengths*ct:2*numGateLengths*(ct+1)]
		chunk = seqs[ct::numRandomizations]
		#Tack on the calibration scalings
		chunk += [[Id(qubit), measBlock], [X(qubit), measBlock]]
		fileNames = compile_to_hardware(chunk, 'RB/RB', suffix='_{0}'.format(ct+1))

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin

def SingleQubitRBT(qubit, seqFileDir, analyzedPulse, showPlot=False):
	"""

	Single qubit randomized benchmarking using atomic Clifford pulses. 

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel) 
	seqFile : file containing sequence strings
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""	
	#Setup a pulse library
	pulseLib = [AC(qubit, cliffNum) for cliffNum in range(24)]
	pulseLib.append(analyzedPulse)
	measBlock = MEAS(qubit)

	seqs = []
	for ct in range(10):
		fileName = 'RBT_Seqs_fast_{0}_F1.txt'.format(ct+1)
		with open(os.path.join(seqFileDir, fileName),'r') as FID:
			fileReader = reader(FID)
			for pulseSeqStr in fileReader:
				seq = []
				for pulseStr in pulseSeqStr:
					seq.append(pulseLib[int(pulseStr)-1])
				seq.append(measBlock)
				seqs.append(seq)

	seqsPerFile = 79
	numFiles = len(seqs)//seqsPerFile

	for ct in range(numFiles):
		chunk = seqs[ct*seqsPerFile:(ct+1)*seqsPerFile]
		#Tack on the calibration scalings
		numCals = 4
		chunk += [[Id(qubit), measBlock]]*numCals + [[X(qubit), measBlock]]*numCals
		fileNames = compile_to_hardware(chunk, 'RBT/RBT', suffix='_{0}'.format(ct+1))

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin








