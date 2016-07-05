from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from ..Cliffords import clifford_seq, clifford_mat, inverse_clifford
from .helpers import create_cal_seqs

import os
from csv import reader
import numpy as np
from functools import reduce

def create_RB_seqs(numQubits, lengths, repeats=32, interleaveGate=None):
	"""
	Create a list of lists of Clifford gates to implement RB.
	"""
	if numQubits == 1:
		cliffGroupSize = 24
	elif numQubits == 2:
		cliffGroupSize = 11520
	else:
		raise Exception("Can only handle one or two qubits.")

	#Create lists of of random integers
	#Subtract one from length for recovery gate
	seqs = []
	for length in lengths:
		seqs += np.random.random_integers(0, cliffGroupSize-1, (repeats, length-1)).tolist()

	#Possibly inject the interleaved gate
	if interleaveGate:
		newSeqs = []
		for seq in seqs:
			newSeqs.append(np.vstack((np.array(seq, dtype=np.int), interleaveGate*np.ones(len(seq), dtype=np.int))).flatten(order='F').tolist())
		seqs = newSeqs
	#Calculate the recovery gate
	for seq in seqs:
		if len(seq) == 1:
			mat = clifford_mat(seq[0], numQubits)
		else:
			mat = reduce(lambda x,y: np.dot(y,x), [clifford_mat(c, numQubits) for c in seq])
		seq.append(inverse_clifford(mat))

	return seqs

def SingleQubitRB(qubit, seqs, showPlot=False):
	"""

	Single qubit randomized benchmarking using 90 and 180 generators.

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	seqs : list of lists of Clifford group integers
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""

	seqsBis = []
	for seq in seqs:
		seqsBis.append(reduce(operator.add, [clifford_seq(c, qubit) for c in seq]))

	#Add the measurement to all sequences
	for seq in seqsBis:
		seq.append(MEAS(qubit))

	#Tack on the calibration sequences
	seqsBis += create_cal_seqs((qubit,), 2)

	fileNames = compile_to_hardware(seqsBis, 'RB/RB')
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)

def TwoQubitRB(q1, q2, seqs, showPlot=False, suffix=""):
	"""

	Two qubit randomized benchmarking using 90 and 180 single qubit generators and ZX90

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	seqs : list of lists of Clifford group integers
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""

	seqsBis = []
	for seq in seqs:
		seqsBis.append(reduce(operator.add, [clifford_seq(c, q1, q2) for c in seq]))

	#Add the measurement to all sequences
	for seq in seqsBis:
		seq.append(MEAS(q1) * MEAS(q2))

	#Tack on the calibration sequences
	seqsBis += create_cal_seqs((q1,q2), 2)

	fileNames = compile_to_hardware(seqsBis, 'RB/RB', suffix=suffix)
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)

def SingleQubitRB_AC(qubit, seqs, showPlot=False):
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
	seqsBis = []
	for seq in seqs:
		seqsBis.append([AC(qubit, c) for c in seq])

	#Add the measurement to all sequences
	for seq in seqsBis:
		seq.append(MEAS(qubit))

	#Tack on the calibration sequences
	seqsBis += create_cal_seqs((qubit,), 2)

	fileNames = compile_to_hardware(seqsBis, 'RB/RB')
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)

def SingleQubitIRB_AC(qubit, seqFile, showPlot=False):
	"""

	Single qubit interleaved randomized benchmarking using atomic Clifford pulses.

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
				seq.append(pulseLib[int(pulseStr)])
			seq.append(measBlock)
			seqs.append(seq)

	#Hack for limited APS waveform memory and break it up into multiple files
	#We've shuffled the sequences so that we loop through each gate length on the inner loop
	numRandomizations = 36
	for ct in range(numRandomizations):
		chunk = seqs[ct::numRandomizations]
		chunk1 = chunk[::2]
		chunk2 = chunk[1::2]
		#Tack on the calibration scalings
		chunk1 += [[Id(qubit), measBlock], [X(qubit), measBlock]]
		fileNames = compile_to_hardware(chunk1, 'RB/RB', suffix='_{0}'.format(2*ct+1))
		chunk2 += [[Id(qubit), measBlock], [X(qubit), measBlock]]
		fileNames = compile_to_hardware(chunk2, 'RB/RB', suffix='_{0}'.format(2*ct+2))

	if showPlot:
		plot_pulse_files(fileNames)

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
		tmpSeqs = []
		with open(os.path.join(seqFileDir, fileName),'r') as FID:
			fileReader = reader(FID)
			for pulseSeqStr in fileReader:
				seq = []
				for pulseStr in pulseSeqStr:
					seq.append(pulseLib[int(pulseStr)-1])
				seq.append(measBlock)
				tmpSeqs.append(seq)
			seqs += tmpSeqs[:12]*12 + tmpSeqs[12:-12] + tmpSeqs[-12:]*12

	seqsPerFile = 100
	numFiles = len(seqs)//seqsPerFile

	for ct in range(numFiles):
		chunk = seqs[ct*seqsPerFile:(ct+1)*seqsPerFile]
		#Tack on the calibration scalings
		numCals = 4
		chunk += [[Id(qubit), measBlock]]*numCals + [[X(qubit), measBlock]]*numCals
		fileNames = compile_to_hardware(chunk, 'RBT/RBT', suffix='_{0}'.format(ct+1))

	if showPlot:
		plot_pulse_files(fileNames)

def SimultaneousRB_AC(qubits, seqs, showPlot=False):
	"""

	Simultaneous randomized benchmarking on multiple qubits using atomic Clifford pulses.

	Parameters
	----------
	qubits : iterable of logical channels to implement seqs on (list or tuple)
	seqs : a tuple of sequences created for each qubit in qubits
	showPlot : whether to plot (boolean)

	Example
	-------
	>>> q1 = QubitFactory('q1')
	>>> q2 = QubitFactory('q2')
	>>> seqs1 = create_RB_seqs(1, [2, 4, 8, 16])
	>>> seqs2 = create_RB_seqs(1, [2, 4, 8, 16])
	>>> SimultaneousRB_AC((q1, q2), (seqs1, seqs2), showPlot=False)
	"""
	seqsBis = []
	for seq in zip(*seqs):
		seqsBis.append([reduce(operator.__mul__, [AC(q,c) for q,c in zip(qubits,
				pulseNums)]) for pulseNums in zip(*seq)])

	#Add the measurement to all sequences
	for seq in seqsBis:
		seq.append(reduce(operator.mul, [MEAS(q) for q in qubits]))

	#Tack on the calibration sequences
	seqsBis += create_cal_seqs((qubits), 2)

	fileNames = compile_to_hardware(seqsBis, 'RB/RB')
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)
