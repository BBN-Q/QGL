from QGL import *

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
	print('Number of sequences: {0}'.format(len(seqs)))

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
	measBlock = MEAS(qubit)

	with open(seqFile,'r') as FID:
		fileReader = reader(FID)
		seqs = []
		for pulseSeqStr in fileReader:
			seq = []
			for pulseStr in pulseSeqStr:
				seq.append(pulseLib[int(pulseStr)])
			seq.append(measBlock)
			seqs.append(seq)

	# #Tack on the calibration scalings
	# seqs += [[Id(qubit), measBlock], [Id(qubit), measBlock], [X(qubit), measBlock], [X(qubit), measBlock]]
	# print('Number of sequences: {0}'.format(len(seqs)))

	# fileNames = compile_to_hardware(seqs, 'RB/RB')
	# print fileNames
	#Hack for limited APS waveform memory and break it up into multiple files
	#We've shuffled the sequences so that we loop through each gate length on the inner loop
	numRandomizations = 32
	numGateLengths = 9
	for ct in range(numRandomizations//2):
		chunk = seqs[2*numGateLengths*ct:2*numGateLengths*(ct+1)]
		#Tack on the calibration scalings
		chunk += [[Id(qubit), measBlock], [Id(qubit), measBlock], [X(qubit), measBlock], [X(qubit), measBlock]]
		print('Number of sequences: {0}'.format(len(chunk)))
		fileNames = compile_to_hardware(chunk, 'RB/RB', suffix='_{0}'.format(ct+1))

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin


