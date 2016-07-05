from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs
from itertools import product
import operator
from ..ControlFlow import *
from functools import reduce

@qfunction
def qreset(qubits, signVec, measDelay, buf):
    # for each qubit, build the set of feedback actions to perform when
    # receiving a zero or one in the comparison register
    FbGates = []
    for ct, q in enumerate(qubits):
        if signVec[ct] == 0:
            FbGates.append([gate(q) for gate in [Id, X]])
        else: # inverted logic
            FbGates.append([gate(q) for gate in [X, Id]])
    FbSeq = [reduce(operator.mul, x) for x in product(*FbGates)]

    # load register
    seq = [Id(qubits[0], measDelay), qwait('CMP'), Id(qubits[0], buf)]
    # create a branch for each possible comparison value
    for ct in range(2**len(qubits)):
        seq += qif(ct, [FbSeq[ct]])

    return seq

def Reset(qubits, measDelay = 1e-6, signVec = None, doubleRound = True, buf = 20e-9, showPlot=False, measChans=None, docals=True, calRepeats=2):
	"""

	Variable amplitude Rabi nutation experiment for an arbitrary number of qubits simultaneously

	Parameters
	----------
	qubits : tuple of logical channels to implement sequence (LogicalChannel)
	measDelay : delay between end of measuerement and LOADCMP
	signVec : conditions for feedback. Tuple of 0 (flip if signal is above threshold) and 1 (flip if below) for each qubit. Default = 0 for all qubits
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
		signVec = (0,)*len(qubits)

	seqs = [prep + [qreset(qubits, signVec, measDelay, buf)] for prep in create_cal_seqs(qubits,1)]
	measBlock = reduce(operator.mul, [MEAS(q) for q in qubits])
	if doubleRound:
		for seq in seqs:
			seq += [measBlock]
			seq.append(qreset(qubits, signVec, measDelay, buf))

	# add final measurement
	for seq in seqs:
		seq += [measBlock, Id(qubits[0], measDelay), qwait('CMP')]

	if docals:
		seqs += create_cal_seqs(qubits, calRepeats, measChans=measChans, waitcmp=True)

	return seqs
