from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from helpers import create_cal_seqs
from itertools import product
import operator
from ..ControlFlow import *

def Reset(qubits, measDelay = 1e-6, signVec = None, doubleRound = True, buf = 30e-9, showPlot=False, measChans=None, docals=True, calRepeats=2):
	"""

	Variable amplitude Rabi nutation experiment for an arbitrary number of qubits simultaneously

	Parameters
	----------
	qubits : tuple of logical channels to implement sequence (LogicalChannel)
	measDelay : delay between end of measuerement and LOADCMP
	signVec : conditions for feedback. List of 0 (flip if signal is above threshold) and 1 (flip if below) for each qubit. Default = 0 for all qubits
	doubleRound : if true, double round of feedback
	showPlot : whether to plot (boolean)
    measChans : tuble of qubits to be measured (LogicalChannel)
    docals, calRepeats: enable calibration sequences, repeated calRepeats times

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	if measChans is None:
		measChans = qubits 

	if signVec == None:
 		signVec = [0]*len(qubits)
        
	states = create_cal_seqs(qubits,1,measChans=measChans)
 	FbSet = [Id, X]
 	FbSet2 = [X, Id]
 	FbGates = []

 	for count in range(len(qubits)):
 		FbGates += [FbSet] if signVec[count]==0 else [FbSet2]    
 	FbSeq = [reduce(operator.mul, [p(q) for p,q in zip(pulseSet, qubits)]) for pulseSet in product(*FbGates)]
	seqs = [state + [MEAS(*measChans), Id(qubits[0],measDelay), qwait('CMP'), Id(qubits[0],buf)] + [branch for b in [qif(fbcount,[FbSeq[count]]) for fbcount in range(len(states))] for branch in b] + [MEAS(*measChans)] for count, state in enumerate(states)] 
    
	if doubleRound:
		seqs = [seq + [Id(qubits[0],measDelay), qwait('CMP'), Id(qubits[0],buf)] + [branch for b in [qif(fbcount,[FbSeq[count]]) for fbcount in range(2**len(qubits))] for branch in b] + [MEAS(*measChans)] for seq in seqs]
	print seqs[0]
	if docals:
		seqs += create_cal_seqs(qubits, calRepeats, measChans=measChans)

	fileNames = compile_to_hardware(seqs, 'Reset/Reset')

	if showPlot:
		plot_pulse_files(fileNames)