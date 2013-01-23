from QGL import *

from csv import reader

def SingleQubitRB(qubit, seqFile, showPlot=False):
	"""

	Single qubit randomized benchmarking. 

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
			tmpSeq = []
			for tmpPulseStr in pulseSeqStr:
				tmpSeq.append(pulseLib[tmpPulseStr])
			tmpSeq.append(measBlock)
			seqs.append(tmpSeq)

	#Tack on the calibration scalings
	seqs += [[Id(qubit), measBlock], [Id(qubit), measBlock], [X(qubit), measBlock], [X(qubit), measBlock]]
	print('Number of sequences: {0}'.format(len(seqs)))

	fileNames = compile_to_hardware(seqs, 'RB/RB')
	print(fileNames)

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin
